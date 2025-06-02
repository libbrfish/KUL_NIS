(ns register
  (:require [util :refer [run-command ensure-dir file-label]]))

(defn apply-transform [input output reference transform interpolation-type]
  (run-command "antsApplyTransforms -d 3 --float 1"
               "--verbose 1"
               "-i" input
               "-o" output
               "-r" reference
               "-n" interpolation-type
               "-t" transform))

(defn rigid-register [source target out-prefix interpolation-type]
  (let [warp-prefix out-prefix
        output-img  (str out-prefix ".nii.gz")]
    (run-command "antsRegistration"
                 "--verbose 1 --dimensionality 3"
                 (str "--output [" warp-prefix "," output-img "]")
                 (str "--interpolation " interpolation-type)
                 "--use-histogram-matching 0"
                 "--winsorize-image-intensities [0.005,0.995]"
                 (str "--initial-moving-transform [" target "," source ",1]")
                 "--transform Rigid[0.1]"
                 "--metric MI[" target "," source ",1,32,Regular,0.25]"
                 "--convergence [1000x500x250x100,1e-6,10]"
                 "--shrink-factors 8x4x2x1"
                 "--smoothing-sigmas 3x2x1x0vox")

    ;; Convert .mat to .txt
    (run-command  "ConvertTransformFile 3" (str warp-prefix "0GenericAffine.mat")
                  (str warp-prefix "0GenericAffine.txt"))))

(defn affine-register-command
  [ants-verbose warp-field output-mri interpolation-type target-mri source-mri]
  (run-command
   "antsRegistration" "--verbose" ants-verbose "--dimensionality" "3"
   "--output" (str "[" warp-field "," output-mri "]")
   "--interpolation" interpolation-type
   "--use-histogram-matching" "0" "--winsorize-image-intensities" "[0.005,0.995]"
   "--initial-moving-transform" (str "[" target-mri "," source-mri ",1]")
   "--transform" "Affine[0.1]"
   "--metric" (str "MI[" target-mri "," source-mri ",1,32,Regular,0.25]")
   "--convergence" "[1000x500x250x100,1e-6,10]"
   "--shrink-factors" "8x4x2x1" "--smoothing-sigmas" "3x2x1x0vox"))

(defn affine-register
  [{:keys [source-mri target-mri ants-verbose interpolation-type label output-dir]
    :or   {ants-verbose       0
           interpolation-type "BSpline"}}]
  (let [registeroutputdir (str output-dir "/antsregister")
        _                 (ensure-dir registeroutputdir)
        label             (or label (file-label source-mri))
        target-mri-label  (file-label target-mri)
        warp-field        (str registeroutputdir "/" label "_reg2_" target-mri-label)
        output-mri        (str output-dir "/" label "_reg2_" target-mri-label ".nii.gz")]
    (affine-register-command ants-verbose
                             warp-field
                             output-mri
                             interpolation-type
                             target-mri
                             source-mri)
    {:warp-field (str warp-field "0GenericAffine.mat")
     :output-mri output-mri}))
