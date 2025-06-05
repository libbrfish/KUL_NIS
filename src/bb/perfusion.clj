(ns perfusion
  (:require
   [babashka.fs :as fs]
   [bet]
   [cheshire.core :as json]
   [clojure.java.io :as io]
   [clojure.tools.logging :as log]
   [register :refer [register]]
   [sc.api]
   [util :refer [find-files run-command with-temp-file]]))

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
  ([dsc-file output]
   (util/ensure-dir (str (babashka.fs/parent (babashka.fs/file output))))
   (with-temp-file [dsc-part ".nii.gz"
                    dsc-avg ".nii.gz"
                    mask ".nii.gz"]
     (get-volumes dsc-file 0 5 dsc-part)
     (register/avg dsc-part dsc-avg)
     (register/motion-correct dsc-file dsc-avg output)
     (denoise output)
     (let [mask-name (str (util/file-label mask) "_bet.nii.gz")
           mask-path (-> mask fs/parent)
           mask-bet  (str mask-path "/" mask-name)]
       (log/info "Generating HDBET mask of DSC series.")
       (bet/hdbet dsc-avg mask true)
       {:output output
        :mask   mask-bet
        :avg    dsc-avg})))
  ([dsc-file output mask-output]
   (util/ensure-dir (str (babashka.fs/parent (babashka.fs/file output))))
   (with-temp-file [dsc-part ".nii.gz"
                    dsc-avg ".nii.gz"
                    mask ".nii.gz"]
     (get-volumes dsc-file 0 5 dsc-part)
     (register/avg dsc-part dsc-avg)
     (register/motion-correct dsc-file dsc-avg output)
     (denoise output)
     (let [mask-name (str (util/file-label mask) "_bet.nii.gz")
           mask-path (-> mask fs/parent)
           mask-bet  (str mask-path "/" mask-name)]
       (log/info "Generating HDBET mask of DSC series.")
       (bet/hdbet dsc-avg mask true)
       (fs/copy mask-bet mask-output)
       {:output output
        :mask   mask-output
        :avg    dsc-avg}))))

(defn calculate-rCBV
  [dsc-reg dir]
  (run-command "src/R/perfusion.R" "-i" (str dsc-reg) "-o" (str dir)))

(defn calculate-rCBV-py
  [dsc-reg mask dir te tr]
  (run-command "src/python/perfusion.py" (str dsc-reg)
               (str mask) (str dir)
               "--te" te
               "--tr" tr))

(defn process-DSC  
  [dir p]
  (let [dsc            (get-native-DSC dir p)
        out-dir        (str (util/kulderivativesdir p dir) "/perfusion")
        t1-bet         (str (util/kulderivativesdir p dir) "/KUL_anat_register/T1w_masked.nii.gz")
        warp-field     (str out-dir "/DSC_reg2_T1")
        transform-file (str warp-field "0GenericAffine.mat")
        rcbv           (str out-dir "/rCBV_corrected.nii.gz")
        rcbv-reg       (str out-dir "/rCBV_reg2_T1.nii.gz")
        rcbf           (str out-dir "/rCBVF.nii.gz")
        rcbf-reg       (str out-dir "/rCBF_reg2_T1.nii.gz")
        mtt            (str out-dir "/MTT.nii.gz")
        mtt-reg        (str out-dir "/MTT_reg2_T1.nii.gz")
        mask-out       (str out-dir "/DSC_mask.nii.gz")]
    (with-temp-file [dsc-preproc ".nii.gz"
                     dsc-reg2-t1 ".nii.gz"
                     dsc-avg-bet ".nii.gz"]
      (let [{:keys [TE TR]} (get-TE-TR dir p)
            dsc-prepared    (prepare-DSC dsc dsc-preproc)
            dsc-mask        (:mask dsc-prepared)
            dsc-avg         (:avg dsc-prepared)]
        #_(calculate-rCBV dsc-preproc out-dir)
        (sc.api/spy (calculate-rCBV-py dsc-preproc dsc-mask out-dir TE TR))
        ;; Copy the mask to the output dir
        (fs/copy dsc-mask mask-out)
        (bet/apply-mask dsc-avg dsc-mask dsc-avg-bet)
        (register {:type               :affine
                   :ants-verbose       0
                   :warp-field         warp-field
                   :output-mri         (str dsc-reg2-t1)
                   :interpolation-type "BSpline"
                   :target-mri         t1-bet
                   :source-mri         dsc-avg-bet})
        (register/apply-transform rcbv rcbv-reg t1-bet transform-file "BSpline")
        (register/apply-transform rcbf rcbf-reg t1-bet transform-file "BSpline")
        (register/apply-transform mtt mtt-reg t1-bet transform-file "BSpline")))))

;; TODO: 1. Masking of vessels on the rCBV images

;; TODO: 2. Normalization of the rCBV maps.

(comment
  (process-DSC "" 1))
