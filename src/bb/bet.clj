(ns bet
  (:require [util :refer [run-command]]))

(defn hdbet
  "mask is saved as out_bet.nii.gz"
  [in out mask?]
  (run-command "hd-bet"
               "-i" in
               "-o" out
               (when mask? "--save_bet_mask")
               "--verbose")
  out)

(defn apply-mask
  [in mask out]
  (run-command
   "mrcalc"
   in mask "-mul" out "-force"))
