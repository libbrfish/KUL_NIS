(ns util
  (:require [clojure.java.io :as io]
            [babashka.process :refer [tokenize process]]
            [clojure.string :as str]
            [clojure.tools.logging :as log]
            [clojure.reflect :as reflect]
            [com.rpl.specter :as s]
            [babashka.fs :as fs]))

;;;; Running commands (shell)

(defn stream-output "Stream output to print."
  [s]
  (with-open [rdr (io/reader s)]
    (binding [*in* rdr]
      (loop []
        (when-let [line (read-line)]
          (log/info line)
          (recur))))))

(defn run-command
  [& commands]
  (let [return-value (->> (str/join " " commands)
                          tokenize
                          (apply process {:continue     true
                                          :pre-start-fn #(apply println "Running: " (:cmd %))
                                          :err          :out}))]
    (-> return-value :out stream-output)
    (let [exit-val (-> return-value deref :exit)]
      (case exit-val
        0 :success
        (str "Returned with exit code " exit-val ". Command: " (str/join " " commands))))))

;;;; Directory commands

(defn ensure-dir [path]
  (let [f (io/file path)]
    (when-not (.exists f)
      (.mkdirs f))))

(defn file-label
  "Get the label of a given file."
  [f]
  (assert f "File-label received nil.")
  (-> f fs/file-name fs/strip-ext fs/strip-ext))

(defn find-file--exclude
  [to-remove v]
  {:pre [(or (string? to-remove)
             (vector? to-remove))]}
  (let [to-remove (if (vector? to-remove) to-remove (vector to-remove))]
    (reduce (fn [v to-remove]
              (filterv #(not (.contains % to-remove)) v))
            v to-remove)))

(defn find-file--include
  [to-include v]
  {:pre [(or (string? to-include)
             (vector? to-include))]}
  (let [to-include (if (vector? to-include) to-include (vector to-include))]
    (reduce (fn [v to-include]
              (filterv #(.contains % to-include) v))
            v to-include)))

(defn find-files
  [dir & {:keys [glob exclude include optional]}]
  (let [file-list (->> (fs/glob dir glob) (mapv str))
        file-list (cond->> file-list
                    (and exclude (seq file-list))
                    (find-file--exclude exclude)
                    
                    (and include (seq file-list))
                    (find-file--include include))]
    (cond
      (not optional) file-list
      :else          (let [file-list-optional (filterv #(.contains % optional) file-list)]
                       (if (seq file-list-optional)
                         file-list-optional
                         file-list)))))

(defn get-bids-dir
  "Find the BIDS dir to work in (recursively from path)."
  [path]
  (-> (fs/match path "glob:**BIDS" {:recursive true})
      first
      fs/absolutize
      str))

(defn kulderivativesdir
  [p dir]
  (str (get-bids-dir dir) "/derivatives/KUL_compute/sub-" p ))

(defn dir-for-participant
  [p dir]
  (let [bids-dir (get-bids-dir dir)
        ;; Is there a biascorrected dir for the subject?
        subj-dir
        (or
         ;; Is there bias corrected data for the subject?
         (seq (fs/match bids-dir (str "glob:**KUL_compute/sub-" p "/KUL_anat_biascorrect")
                        {:recursive true}))
         ;; Is there a subject data dir for the subject?
         (seq (fs/match bids-dir (str "glob:**sub-" p) {:recursive true}))
         (throw (Exception. "Subject data dir not found.")))]
    (-> subj-dir first str)))

(defn images-for-participant
  [p dir]
  (let [sub-dir  (dir-for-participant p dir)
        bids-dir (str  (get-bids-dir dir) "/sub-" p "/")]
    (println sub-dir)
    {:T1w    (find-files sub-dir
                         :glob "**T1w*.nii.gz"
                         :exclude ["gadolinium"
                                   "cT1w"
                                   "_reg2_"]
                         :optional "bc")
     :cT1w   (find-files sub-dir
                         :glob "**cT1w*.nii.gz"
                         :optional "bc")
     :FLAIR  (find-files sub-dir
                         :glob "**FLAIR*.nii.gz"
                         :optional "bc")
     :FGATIR (find-files sub-dir
                         :glob "**FGATIR*.nii.gz"
                         :optional "bc")
     :T2w    (find-files sub-dir
                         :glob "**T2w*.nii.gz"
                         :optional "bc")
     :SWI    (find-files bids-dir
                         :glob "**run-01_SWI*.nii.gz"
                         :optional "bc")
     :SWIp   (find-files bids-dir
                         :glob "**run-02_SWI*.nii.gz"
                         :optional "bc")}))

;;;; Files
(defn- with-temp-file--build-vector [v]
  (loop [v      (->> v reverse (into '()))
         result []]
    
    (if (seq v)
      (recur (-> v pop pop)
             (conj result (peek v) `(fs/create-temp-file {:suffix ~(-> v pop peek)})))
      result)))

(defmacro with-temp-file
  "Execute code with temporary file(s).
  v is a vector with name of var and suffix of temp file."
  {:clj-kondo/lint-as 'clojure.core/let}
  [v & commands]
  `(let ~(with-temp-file--build-vector v)
     (let [return# (do ~@commands)]
       ~@(mapv #(list `fs/delete-on-exit (first %))
               (partition 2 v))
       return#)))


;;;; Reflect
(defn reflect-item
  "What functions can we use on item?"
  [item]
  (s/select [s/ALL :name] (-> item reflect/reflect :members)))
