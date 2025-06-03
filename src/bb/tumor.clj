(ns tumor
  (:require [util :refer [run-command]]
            [register :refer [affine-register-command apply-transform]]
            [babashka.fs :as fs]))

(defn preproc-images
  [p dir]
  (let [ preproc-dir (str  (util/kulderivativesdir p dir) "/KUL_anat_register")]
    {:T1w   (str preproc-dir "/T1w_masked.nii.gz")
     :cT1w  (str preproc-dir "/cT1w_bc_reg2_T1w_masked.nii.gz")
     :FLAIR (str preproc-dir "/FLAIR_bc_reg2_T1w_masked.nii.gz")
     :T2w   (str preproc-dir "/T2w_bc_reg2_T1w_masked.nii.gz")}))

(defn process-tumor
  [p dir]
  (let [tumor-dir                    (str  (util/kulderivativesdir p dir) "/tumor")
        {:keys [T1w cT1w FLAIR T2w]} (preproc-images p dir)
        mni-image                    "/opt/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"
        warp-field                   (fs/create-temp-file)
        interpolation-type           "BSpline"]
    (affine-register-command 0
                             warp-field
                             (str tumor-dir "/T1.nii.gz")
                             "BSpline"
                             mni-image
                             T1w)
    (let [transform (str warp-field "0GenericAffine.mat")
          tumor-mni (str tumor-dir "/tumor_MNI.nii.gz")
          tumor-sub (str tumor-dir "/tumor_subj-space.nii.gz")]
      (apply-transform T2w (str tumor-dir "/T2.nii.gz") mni-image transform interpolation-type)
      (apply-transform cT1w (str tumor-dir "/CT1.nii.gz") mni-image transform interpolation-type)
      (apply-transform FLAIR (str tumor-dir "/FLAIR.nii.gz") mni-image transform interpolation-type)
      (run-command "hd_glio_predict" "-t1" (str tumor-dir "/T1.nii.gz")
                   "-t1c" (str tumor-dir "/CT1.nii.gz")
                   "-t2" (str tumor-dir "/T2.nii.gz")
                   "-flair" (str tumor-dir "/FLAIR.nii.gz")
                   "-o" tumor-mni)
      ;; Reverse transvorm 
      (apply-transform tumor-mni tumor-sub T1w transform "NearestNeighbor" true))))
