(ns perfusion
  (:require [util :refer [find-files run-command
                          with-temp-file]]
            [register :refer []]))

(defn get-native-DSC
  [dir p]
  (first (find-files dir
                     :glob (str "**sub-" p "**dsc.nii.gz"))))

(defn get-volumes
  [img-4d start end o]
  (run-command
   "mrconvert" (str img-4d)  (str o)
   "-coord 3" (str start ":" end)
   "-force"))


(defn prepare-DSC
  [dsc-file output]
  (util/ensure-dir (str (babashka.fs/parent (babashka.fs/file output))))
  (with-temp-file [dsc-part ".nii.gz"
                   dsc-avg ".nii.gz"]
    (get-volumes dsc-file 0 5 dsc-part)
    (register/avg dsc-part dsc-avg)
    (register/motion-correct dsc-file dsc-avg output)))
