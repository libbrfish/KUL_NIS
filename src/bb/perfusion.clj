(ns perfusion
  (:require [util :refer [find-files run-command
                          with-temp-file]]
            [register :refer [register]]
            [cheshire.core :as json]
            [clojure.java.io :as io]))

(defn get-native-DSC
  [dir p]
  (first (find-files dir :glob (str "**sub-" p "**dsc.nii.gz"))))

(defn get-native-DSC-json
  [dir p]
  (first (find-files dir :glob (str "**sub-" p "**dsc.json"))))

(defn get-TE--
  [m]
  (-> m :EchoTime (* 1000)))

(defn get-TR--
  [m]
  (-> m :RepetitionTime (* 1000)))

(defn get-TE-TR
  [dir p]
  (let [dsc-json 
        (->> (get-native-DSC-json dir p)
             io/reader
             json/parse-stream 
             (mapv (fn [[k v]] [(keyword k) v]))
             (into {}))]
    {:TE (get-TE-- dsc-json)
     :TR (get-TR-- dsc-json)}))

(defn get-volumes
  [img-4d start end o]
  (run-command
   "mrconvert" (str img-4d)  (str o)
   "-coord 3" (str start ":" end)
   "-force"))

(defn denoise
  [dsc]
  (run-command "dwidenoise" dsc dsc "-force"))

(defn prepare-DSC
  "Motion correction and denoising of DSC image data."
  [dsc-file output]
  (util/ensure-dir (str (babashka.fs/parent (babashka.fs/file output))))
  (with-temp-file [dsc-part ".nii.gz"
                   dsc-avg ".nii.gz"]
    (get-volumes dsc-file 0 5 dsc-part)
    (register/avg dsc-part dsc-avg)
    (register/motion-correct dsc-file dsc-avg output)
    (denoise output)
    output))

(defn calculate-rCBV
  [dsc-reg dir]
  (run-command "src/R/perfusion.R" "-i" (str dsc-reg) "-o" (str dir)))

(defn process-DSC
  [dir p]
  (let [dsc            (get-native-DSC dir p)
        out-dir        (str (util/kulderivativesdir p dir) "/perfusion")
        t1-bet         (str (util/kulderivativesdir p dir) "/tumor/T1.nii.gz")
        dsc-avg-bet    (str out-dir "/R_avg_bet.nii.gz")
        warp-field     (str out-dir "/DSC_reg2_T1")
        transform-file (str warp-field "0GenericAffine.mat")
        rcbv           (str out-dir "/rCBV_corrected.nii.gz")
        rcbv-reg       (str out-dir "/rCBV_reg2_T1.nii.gz")]
    (with-temp-file [dsc-preproc ".nii.gz"
                     dsc-reg2-t1 ".nii.gz"]
      (prepare-DSC dsc dsc-preproc)
      (calculate-rCBV dsc-preproc out-dir)
      (register {:type               :affine
                 :ants-verbose       0
                 :warp-field         warp-field
                 :output-mri         (str dsc-reg2-t1)
                 :interpolation-type "BSpline"
                 :target-mri         t1-bet
                 :source-mri         dsc-avg-bet})
      (register/apply-transform rcbv rcbv-reg t1-bet transform-file "BSpline"))))

;; TODO: 1. Masking of vessels on the rCBV images

;; TODO: 2. Normalization of the rCBV maps.
