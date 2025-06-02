(ns perfusion-spec
  (:require [speclj.core :refer :all]
            [babashka.fs :as fs]
            [perfusion :refer :all]
            [clojure.reflect :as reflect]
            [util :refer [reflect-item]]))

(describe "get-native-DSC"
          (it "Can find the absolute file path of the DSC nii file."
              (should=
               "/home/devos/Documents/projects/DSC_pipeline/resources/BIDS/sub-1/perf/sub-1_dsc.nii.gz"
               (get-native-DSC "/home/devos/Documents/projects/DSC_pipeline/" 1))))

(describe "prepare-DSC"
          (with-stubs)
          (it "Performs motion correction on native DSC images"
              (with-redefs
                [util/run-command (fn [& x]
                                    (println x)
                                    (clojure.string/join " " x))]
                (let [result (prepare-DSC
                              "DSC_in.nii.gz"
                              "DSC_mc.nii.gz")]
                  (should-contain "antsMotionCorr --dimensionality 3 --n-images 10 --metric"
                                  result)
                  (should-contain ".nii.gz,DSC_in.nii.gz,1,32,Regular,0.5] --useFixedReferenceImage 1 --transform Rigid[0.2] --iterations 4x2x1 --smoothingSigmas 4x2x0 --shrinkFactors 8x4x1 --output [corr_,DSC_mc.nii.gz"
                                  result)))))

(describe "calculate-rCBV"
          (it "Can calculate the rCBV maps from the registered DSC images."
              (should= :succes
                       (calculate-rCBV "DSC-reg.nii.gz" "DSC-rCBV.nii.gz"))))

(run-specs)
