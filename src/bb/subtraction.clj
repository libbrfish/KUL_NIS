(ns subtraction
  (:require [util]
            [bet]
            [clojure.string :as str]
            [babashka.fs :as fs]
            [clojure.tools.logging :as log]))

(defn mrstats
  "Possible output: mean, median, std, std_rv, min, max, count"
  [image output mask]
  (-> (babashka.process/sh "mrstats" "-mask" mask "-output" output image)
      :out (str/split #" ") first Float/parseFloat))

(defn gen-p95-mask
  [image mask output]
  (let [max (mrstats image "max" mask)
        min (mrstats image "min" mask)
        p95 (-> (- max min)
                (/ 100)
                (* 35)
                (+ min))]
    (util/run-command "mrcalc" image p95 "-lt" output "-force")))

(defn gen-histogram-mask
  "Generate a mask to use for the histogram matching.
This mask should be comprised of:
- tumor and peritumoral edema with some dilatation.
- BET mask.
- Masking of the >p95 structures on the BET post contrast images."
  [bet-mask tumor-mask cT1 output]
  (util/with-temp-file [tumor-mask-dilated ".nii.gz"
                        bet-tumor-mask ".nii.gz"
                        cT1-masked ".nii.gz"
                        p95-mask ".nii.gz"]
    ;; Dilatate the tumor mask
    (util/run-command "maskfilter" "-force"
                      "-npass 5"
                      tumor-mask
                      "dilate"
                      tumor-mask-dilated)
    ;; Calculate the combination of the bet-mask and tumor-mask-dilated
    (util/run-command "mrcalc" tumor-mask-dilated 0 "-eq"
                      bet-mask "-min"
                      bet-tumor-mask "-force")
    
    ;; Apply this mask to the cT1 images
    (util/run-command "mrcalc" cT1 bet-tumor-mask "-mult" cT1-masked "-force")

    ;; Mask the top 5% intensity values of the cT1 images.
    (log/info "Running gen-p95-mask.")
    (gen-p95-mask cT1 bet-tumor-mask p95-mask)
    
    ;; Combine bet-mask, tumor-mask and p95-mask
    (util/run-command "mrcalc" p95-mask bet-tumor-mask "-min" output "-force")
    output))

(defn get-bet-mask
  "Get the BET mask, which should be located in the KUL_anat_register dir."
  [p dir]
  (let [KUL_anat_register-dir (str (util/kulderivativesdir p dir) "/KUL_anat_register")
        bet-mask              (str KUL_anat_register-dir "/T1w_masked_bet.nii.gz")]
    bet-mask))

(defn get-tumor-mask
  "Get the tumor mask, which should be located in the tumor dir."
  [p dir]
  (let [tumor-dir  (str (util/kulderivativesdir p dir) "/tumor")
        tumor-mask (str tumor-dir "/tumor_subj-space.nii.gz")]
    tumor-mask))

(defn match-histogram
  "Match the histogram of source image to target image to produce output image."
  [source target output bet-mask tumor-mask]
  (util/with-temp-file [mask ".nii.gz"]
    (let [cT1  target
          T1   source
          mask (str mask)]
      (gen-histogram-mask bet-mask tumor-mask cT1 mask)
      (util/run-command "mrhistmatch nonlinear" T1 cT1 output
                        "-bins 5"
                        "-mask_input" mask
                        "-mask_target" mask
                        "-force")
      mask)))

(defn process
  [p dir]
  (let [kul-dir     (util/kulderivativesdir p dir)
        output-dir  (str kul-dir "/subtraction")
        _           (util/ensure-dir output-dir)
        source      (str kul-dir "/KUL_anat_register/T1w_masked.nii.gz")
        target      (str kul-dir "/KUL_anat_register/cT1w_bc_reg2_T1w_masked.nii.gz")
        output      (str output-dir "/T1w_hm2_cT1w.nii.gz")
        subtraction (str output-dir "/subtraction.nii.gz")
        bet-mask    (get-bet-mask p dir)
        tumor-mask  (get-tumor-mask p dir)
        total-mask  (match-histogram source target output bet-mask tumor-mask)]
    (fs/copy target (str output-dir "/cT1w.nii.gz"))
    (fs/move total-mask (str output-dir "/mask.nii.gz"))
    (util/run-command "mrcalc" target output "-sub" subtraction "-force")))
