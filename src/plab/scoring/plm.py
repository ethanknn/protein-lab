"""Protein language model scoring with ESM-2.

What this measures
------------------
ESM-2 was trained on roughly 250 million natural protein sequences. Asking it
how surprised it is by a candidate is a genuine, independent check on whether
the sequence looks like something evolution could have produced -- which
correlates with being foldable, expressible and stable.

The metric is masked-marginal pseudo-perplexity: mask each position in turn,
ask the model for the probability of the residue that is actually there, and
take the exponential of the mean negative log-likelihood. Lower is more
protein-like. This is the standard formulation from the ESM papers.

What it does not measure: binding. ESM-2 has never seen our target. A low
perplexity means "plausible protein", not "good binder". It is weighted in
the composite accordingly.
"""

from __future__ import annotations

import math
import os
from functools import lru_cache

from .base import ScoreComponent, Scorer, ramp_score

# Default to the 650M model: it fits comfortably in 8 GB of VRAM and is the
# usual accuracy/speed compromise. Override with PLAB_ESM_MODEL.
DEFAULT_MODEL = os.environ.get("PLAB_ESM_MODEL", "esm2_t33_650M_UR50D")

# Reference points for turning perplexity into a 0..1 score. Natural globular
# domains typically land in the 5-12 range for this model; shuffled sequences
# of the same composition sit far higher. Calibrated per-run against the
# reference binder, so these are only fallbacks.
PPL_GOOD = 8.0
PPL_POOR = 22.0


class _EsmRunner:
    """Lazily loads ESM-2 and caches scores. One instance per process."""

    def __init__(self, model_name: str = DEFAULT_MODEL) -> None:
        self.model_name = model_name
        self._model = None
        self._alphabet = None
        self._converter = None
        self._device = None
        self._load_error: str | None = None
        self._cache: dict[str, float] = {}

    # --- loading ---------------------------------------------------------

    def load(self) -> bool:
        if self._model is not None:
            return True
        if self._load_error is not None:
            return False
        try:
            import torch
            import esm
        except ImportError as exc:
            self._load_error = f"torch/fair-esm not importable: {exc}"
            return False
        try:
            loader = getattr(esm.pretrained, self.model_name)
            model, alphabet = loader()
            model.eval()
            device = "cuda" if torch.cuda.is_available() else "cpu"
            model = model.to(device)
            self._model = model
            self._alphabet = alphabet
            self._converter = alphabet.get_batch_converter()
            self._device = device
            self._torch = torch
        except Exception as exc:  # download failure, OOM, etc.
            self._load_error = f"{type(exc).__name__}: {exc}"
            return False
        return True

    @property
    def status(self) -> str:
        if self._model is not None:
            return f"{self.model_name} on {self._device}"
        return self._load_error or "not loaded"

    @property
    def device(self) -> str:
        return self._device or "unavailable"

    # --- scoring ---------------------------------------------------------

    def pseudo_perplexity(self, sequence: str, batch_size: int = 16) -> float | None:
        """Masked-marginal pseudo-perplexity. Returns None if unavailable."""
        if sequence in self._cache:
            return self._cache[sequence]
        if not self.load():
            return None

        torch = self._torch
        model, alphabet, converter = self._model, self._alphabet, self._converter
        mask_idx = alphabet.mask_idx

        _, _, tokens = converter([("q", sequence)])
        tokens = tokens.to(self._device)
        n = len(sequence)
        # Token layout is [cls] + residues + [eos], so residue i is at i+1.
        positions = list(range(1, n + 1))

        total_nll = 0.0
        try:
            with torch.no_grad():
                for start in range(0, len(positions), batch_size):
                    chunk = positions[start:start + batch_size]
                    batch = tokens.repeat(len(chunk), 1)
                    true_ids = []
                    for row, pos in enumerate(chunk):
                        true_ids.append(int(batch[row, pos].item()))
                        batch[row, pos] = mask_idx
                    out = model(batch, repr_layers=[], return_contacts=False)
                    logits = out["logits"]
                    logprobs = torch.log_softmax(logits, dim=-1)
                    for row, pos in enumerate(chunk):
                        total_nll -= float(logprobs[row, pos, true_ids[row]].item())
        except Exception:
            # Out of memory or a backend error: report unavailable rather than
            # silently returning a misleading number.
            return None

        ppl = math.exp(total_nll / max(n, 1))
        self._cache[sequence] = ppl
        return ppl


@lru_cache(maxsize=1)
def get_runner(model_name: str = DEFAULT_MODEL) -> _EsmRunner:
    return _EsmRunner(model_name)


class PLMPlausibilityScorer(Scorer):
    """Is this sequence something evolution plausibly could have made?"""

    name = "plm_plausibility"
    weight = 1.6

    def __init__(self, model_name: str = DEFAULT_MODEL) -> None:
        self.runner = get_runner(model_name)

    def available(self) -> bool:
        return self.runner.load()

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        ppl = self.runner.pseudo_perplexity(sequence)

        if ppl is None:
            return ScoreComponent(
                name=self.name,
                raw=float("nan"),
                normalised=0.5,
                weight=0.0,  # excluded from the composite rather than guessed
                units="pseudo-perplexity",
                detail=(
                    f"ESM-2 unavailable ({self.runner.status}); this component "
                    f"was excluded from the composite rather than estimated."
                ),
                caveat="Component disabled.",
            )

        # Anchor the scale to THIS scaffold's own family, not to the reference
        # binder.
        #
        # Measured on esm2_t33_650M: the 50C1 antibody scFv sits at 2.91 while
        # the wild-type Z domain -- a real, natural protein -- sits at 10.08.
        # ESM-2 has seen enormous numbers of antibodies and comparatively few
        # three-helix bundles, so anchoring an affibody's scale to an antibody
        # compares across protein families rather than measuring design
        # quality. It rated the natural parent protein at 0%, which is a
        # self-evident failure of the scale, and because every candidate then
        # saturated at 0 the component carried weight while supplying no
        # gradient at all -- worse than omitting it.
        #
        # The fair question is: closer to the natural parent, or to a random
        # library member? Both anchors are measured at runtime per scaffold.
        good = context.get("plm_good")
        poor = context.get("plm_poor")
        source = "scaffold-anchored"

        if not good or not poor or poor <= good:
            ref_ppl = context.get("reference_perplexity")
            if ref_ppl:
                # Only sane when the candidate IS the reference's own family.
                good, poor = ref_ppl * 1.05, ref_ppl * 2.6
                source = "reference-anchored"
            else:
                good, poor = PPL_GOOD, PPL_POOR
                source = "fallback scale"

        norm = ramp_score(ppl, poor, good)  # inverted: lower ppl scores higher

        return ScoreComponent(
            name=self.name,
            raw=ppl,
            normalised=norm,
            weight=self.weight,
            units="pseudo-perplexity",
            detail=(
                f"ESM-2 pseudo-perplexity {ppl:.2f} ({self.runner.status}). "
                f"{source}: {good:.2f} = the scaffold's natural parent "
                f"sequence, {poor:.2f} = a random library member on the same "
                f"scaffold. Lower means the model finds the sequence more "
                f"protein-like."
            ),
            caveat=(
                "Measures plausibility as a protein, NOT binding -- ESM-2 has "
                "never seen this target. It also rewards predictability: a "
                "poly-arginine face scores 10.45 here, near the natural "
                "protein's 10.08, despite being a useless binder. Only the "
                "complexity and recognition components catch that, so this "
                "score is never meaningful on its own."
            ),
        )
