(ns register
  (:require [util :refer [run-command ensure-dir file-label]]))

(defmulti register :type)

(defn apply-transform [input output reference transform interpolation-type]
  (run-command "antsApplyTransforms -d 3 --float 1"
               "--verbose 1"
               "-i" input
               "-o" output
               "-r" reference
               "-n" interpolation-type
               "-t" transform))

(defn rigid-register-command
  [ants-verbose warp-field output-mri interpolation-type target-mri source-mri]
  (run-command
   "antsRegistration" "--verbose" ants-verbose "--dimensionality" "3"
   "--output" (str "[" warp-field "," output-mri "]")
   "--interpolation" interpolation-type
   "--use-histogram-matching 0"
   "--winsorize-image-intensities [0.005,0.995]"
   "--initial-moving-transform" (str "[" target-mri "," source-mri ",1]")
   "--transform Rigid[0.1]"
   "--metric" (str "MI[" target-mri "," source-mri ",1,32,Regular,0.25]")
   "--convergence" "[1000x500x250x100,1e-6,10]"
   "--shrink-factors 8x4x2x1"
   "--smoothing-sigmas 3x2x1x0vox"))

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
   "--shrink-factors" "8x4x2x1"
   "--smoothing-sigmas" "3x2x1x0vox"))

(defmethod register :affine
  [{:keys [ants-verbose warp-field output-mri interpolation-type target-mri source-mri]}]
  (affine-register-command ants-verbose warp-field output-mri interpolation-type target-mri source-mri))

(defmethod register :rigid
  [{:keys [ants-verbose warp-field output-mri interpolation-type target-mri source-mri]}]
  (rigid-register-command ants-verbose warp-field output-mri interpolation-type target-mri source-mri))

(defn bids-register
  [{:keys [type source-mri target-mri ants-verbose interpolation-type label output-dir]
    :or   {ants-verbose       0
           interpolation-type "BSpline"}}]
  (let [registeroutputdir (str output-dir "/antsregister")
        _                 (ensure-dir registeroutputdir)
        label             (or label (file-label source-mri))
        target-mri-label  (file-label target-mri)
        warp-field        (str registeroutputdir "/" label "_reg2_" target-mri-label)
        output-mri        (str output-dir "/" label "_reg2_" target-mri-label ".nii.gz")]
    (register {:type               type
               :ants-verbose       ants-verbose
               :warp-field         warp-field
               :output-mri         output-mri
               :interpolation-type interpolation-type
               :target-mri         target-mri
               :source-mri         source-mri})
    {:warp-field (str warp-field "0GenericAffine.mat")
     :output-mri output-mri}))



(defn motion-correct
  "Motion correct a 4d MRI sequence to output o (output is a suffix).
  This will also generate an average image."
  [i o-avg o-corr]
  (run-command
   "antsMotionCorr"
   "--dimensionality 3"
   "--n-images 10"
   "--metric" (str "MI[" o-avg "," i ",1,32,Regular,0.5]")
   "--useFixedReferenceImage 1"
   "--transform Rigid[0.2]"
   "--iterations 4x2x1"
   "--smoothingSigmas 4x2x0"
   "--shrinkFactors 8x4x1"
   "--output" (str "[corr_," o-corr "," o-avg "]")
   "--interpolation BSpline[3] --verbose 1"))

(defn extract-time-points
  "Extract time points a to b (starting with 0) from series i, output to o."
  [i a b o]
  (run-command
   "mrconvert" i "-coord 3" (str a ":1:" b) o "-force")
  o)

(defn avg
  "Calculate the average image o of input i."
  [i o]
  (run-command
   "mrmath" (str i) "mean" (str o) "-axis 3 -force"))
