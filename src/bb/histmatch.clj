(ns histmatch
  (:require
   [babashka.process :as sh]
   [cheshire.core :as json]
   [clojure.string :as str]
   [com.rpl.specter :as s]
   [util]))

(defn get-DSC-data
  [dir p]
  {:rCBV-corrected
   (first (util/find-files dir :glob (str "**sub-" p "**rCBV_corrected.nii.gz")))

   :DSC-mask
   (first (util/find-files dir :glob (str "**sub-" p "**DSC_mask.nii.gz")))})

(comment
  (get-DSC-data "/home/devos/Documents/projects/DSC_pipeline/resources/" 2))

(defn get-percentiles
  [dir p]
  (let [{:keys [rCBV-corrected DSC-mask]} (get-DSC-data dir p)
        percentile-map
        (-> (sh/sh "src/python/percentile.py" rCBV-corrected "-m" DSC-mask)
            :out
            json/parse-string)

        percentile-map
        (->> percentile-map
             (s/transform [s/MAP-KEYS ] keyword)
             (into (sorted-map)))]
    percentile-map))

(defn histmatch
  ([in mask in-hist ref-hist out]
   
   (let [in-hist  (str "\"" (str/join "," (s/select [s/MAP-VALS] in-hist)) "\"")
         ref-hist (str "\"" (str/join "," (s/select [s/MAP-VALS] ref-hist)) "\"")]
     (sh/sh "src/R/histmatch.R"
            "--input" in
            "--output" out
            "--src" in-hist
            "--ref" ref-hist
            "--mask" mask)))
  ([in in-hist ref-hist out]
   (let [in-hist  (str/join "," (s/select [s/MAP-VALS] in-hist))
         ref-hist (str/join "," (s/select [s/MAP-VALS] ref-hist))]
     (util/run-command "src/R/histmatch.R"
                       "--input" in
                       "--output" out
                       "--src" in-hist
                       "--ref" ref-hist))))


;; TODO: Nog om te zetten in code
(comment
  (let [dir                  "/home/devos/Documents/projects/DSC_pipeline/resources/"
        {in      :rCBV-corrected
         in-mask :DSC-mask}  (get-DSC-data dir 1)
        {ref      :rCBV-corrected
         ref-mask :DSC-mask} (get-DSC-data dir 2)
        in-hist              (get-percentiles dir 1)
        ref-hist             (get-percentiles dir 2)
        out                  "/home/devos/Documents/projects/DSC_pipeline/resources/BIDS/derivatives/KUL_compute/sub-1/histogram/rCBV_hm2sub2.nii.gz"]
    (println "in:" in)
    (println "in-mask:" in-mask)
    (println "ref:" ref)
    (println "ref-mask:" ref-mask)
    (println "in-hist:" in-hist)
    (println "ref-hist:" ref-hist)
    (histmatch in in-mask in-hist ref-hist out)))

