(ns spec.bb.KUL-anat-register-spec
  (:require [speclj.core :refer :all]
            [babashka.fs :as fs]
            [KUL-anat-register :refer :all]))

(describe "find-file--remove"
          (it "removes multiple matches."
              (should=
               ["ddd"]
               (find-file--exclude ["aa" "bb"]
                                   ["aaa" "bbb" "aabbcc" "ddd"]))))

(describe "Find-files function"
          (with-stubs)
          (it "Can use glob to detect simple files"
              (with-redefs  [fs/glob (stub :glob {:return ["T1w.nii.gz"
                                                           "gadoliniumT1w.nii.gz"
                                                           "T1w-bc.nii.gz"]})]
                (should= ["T1w.nii.gz" "T1w-bc.nii.gz"] (find-files "/home/devos"
                                                                    :glob "**T1w.nii.gz"
                                                                    :exclude "gadolinium"))
                (should= ["T1w.nii.gz" "T1w-bc.nii.gz"] (find-files "/home/devos"
                                                                    :glob "**T1w.nii.gz"
                                                                    :exclude ["gadolinium"
                                                                              "_reg2_"]))
                (should= ["gadoliniumT1w.nii.gz"] (find-files "/home/devos"
                                                              :glob "**T1w.nii.gz"
                                                              :include ["gadolinium" "nii"]))
                (should= ["T1w-bc.nii.gz"] (find-files "/home/devos"
                                                       :glob "**T1w.nii.gz"
                                                       :exclude "gadolinium"
                                                       :optional "bc")))))

(describe "get-bids-dir"
          (it "Can find the BIDS directory recursively"
              (should= "/home/devos/Documents/projects/DSC_pipeline/resources/BIDS" (get-bids-dir "/home/devos/Documents/projects/DSC_pipeline/resources/"))))

(describe "get-bias-corrected-dir"
          (it "Can find the biascorrect dir of the participant"
              (should=
               "/home/devos/Documents/projects/DSC_pipeline/resources/BIDS/derivatives/KUL_compute/sub-1/KUL_anat_biascorrect"
               (get-bias-corrected-dir 1))))
