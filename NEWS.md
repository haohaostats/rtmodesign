# rtmodesign 1.0.0

* Added the ordinary-user `rtmo_design()` workflow for patient-level pilot
  data from active-controlled dose-response studies.
* Added five built-in continuous dose-response models.
* Added automatic residual moments through order six and adaptive moment
  regularization.
* Added approximate RTMO optimization, sensitivity-function certification,
  and exact floor-ceiling allocation.
* Aligned the internal support search with the manuscript algorithm by using
  the three-point start for the manuscript model, a 10,001-point scan,
  continuous support refinement, and the 1e-8 equivalence threshold.
* Added standard print, summary, plot, and data-frame methods.
