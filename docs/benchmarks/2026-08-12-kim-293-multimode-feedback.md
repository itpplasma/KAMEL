# KIM #293: multi-mode periodic feedback

The periodic adapter now retains signed mode identity (`m`, `n`), solve status,
resonant radius, current-normalization metadata, and compact-transition state
per configured mode. `kim_run_for_all_modes` resets all per-mode arrays before
each batch, refreshes the injected QL-Balance profiles before solving, and
stores each mode independently. `get_dql` then adds the already embedded
per-mode tensors incoherently on the global radial grid.

No absolute-value or mode-index rewrite is performed. A failed solve terminates
the batch with its KIM status; a successful periodic mode carries positive
resonance metadata. Nonperiodic modes retain the initialized zero resonance.
Repeated batches therefore cannot reuse fields, grids, currents, or tensors
from an earlier mode list.
