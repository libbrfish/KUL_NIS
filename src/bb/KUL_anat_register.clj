(ns KUL-anat-register
  (:require [babashka.cli :as cli]
            [babashka.process :refer [process shell sh]]
            [babashka.fs :as fs]
            [clojure.java.io :as io]
            [clojure.tools.logging :as log]
            [clojure.string :as str]))

(def cli-spec
  {:p {:desc "Participant name"}
   :t {:desc "Target image"}
   :s {:desc "Source image"}
   :d {:desc "Output directory"}
   :i {:desc "Interpolation type (1=BSpline, 2=NearestNeighbor, 3=Linear)"}
   :o {:desc "Apply transformation to other images (string or quoted list)"}
   :m {:desc "Mask mode"}
   :r {:desc "Registration type (1=rigid, 2=affine, 3=nonrigid)"}
   :n {:desc "Number of CPUs" :default "15"}
   :v {:desc "Verbosity (0=silent, 1=normal, 2=verbose)" :default "1"}
   :c {:desc "Use bias-corrected input?" :flag true}
   :w {:desc "Use warp2MNI?" :flag true}})

(def args (cli/parse-args *command-line-args* {:spec cli-spec}))



;;;; Utility functions

(defn usage []
  (println "Usage: anat_register.clj -p <participant> | -s <source> -t <target> [options]")
  (println "\nRequired:")
  (println "  -p: participant name (for full BIDS registration)")
  (println "  OR")
  (println "  -s: source image")
  (println "  -t: target image")
  (System/exit 1))

(defn ensure-dir [path]
  (let [f (io/file path)]
    (when-not (.exists f)
      (.mkdirs f))))

(def interpolation-map {"1" "BSpline" "2" "NearestNeighbor" "3" "Linear"})

(def interpolation-type
  (get interpolation-map (str (:i args "1")) "BSpline"))

(defn run [cmd]
  (log/info "Running:" cmd)
  (shell {:out :inherit :err :inherit} cmd))

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

;;;; Directory fns
(defn get-bids-dir
  "Find the BIDS dir to work in (recursively from path)."
  [path]
  (-> (fs/match path "glob:**BIDS" {:recursive true})
      first
      str))

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

(defn file-label
  "Get the label of a given file."
  [f]
  (-> f fs/file-name fs/strip-ext))

(def participant (:p args))
(def kul-log-dir (str "KUL_LOG/KUL_anat_register/sub-" participant))
(ensure-dir kul-log-dir)
(assert (get-bids-dir "."))
(def kulderivativesdir (str (get-bids-dir ".") "derivatives/KUL_compute/sub-" participant "/KUL_anat_register"))
(def registeroutputdir (str kulderivativesdir "/antsregister"))
(def verbose-level (or (:v args) 0))


;;;; ANTs Commands
(defn apply-transform [input output reference transform]
  (run (str/join " " ["antsApplyTransforms -d 3 --float 1"
                      "--verbose 1"
                      "-i" input
                      "-o" output
                      "-r" reference
                      "-n" interpolation-type
                      "-t" transform])))

(defn rigid-register [source target out-prefix]
  (let [warp-prefix out-prefix
        output-img  (str out-prefix ".nii.gz")]
    (run (str/join " " ["antsRegistration"
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
                        "--smoothing-sigmas 3x2x1x0vox"]))

    ;; Convert .mat to .txt
    (run (str "ConvertTransformFile 3 " warp-prefix "0GenericAffine.mat "
              warp-prefix "0GenericAffine.txt"))))

(defn affine-register
  [{:keys [source-mri target-mri ants-verbose interpolation-type]
    :or   {ants-verbose       0
           interpolation-type "BSpline"}}]
  (let [ source-mri-label (file-label source-mri)
        target-mri-label  (file-label target-mri)
        warp-field        (str registeroutputdir "/" source-mri-label "_reg2_" target-mri-label)
        output-mri        (str kulderivativesdir "/" source-mri-label "_reg2_" target-mri-label ".nii.gz")]
    (shell "antsRegistration" "--verbose" ants-verbose "--dimensionality" "3"
           "--output" (str "[" warp-field "," output-mri "]")
           "--interpolation" interpolation-type
           "--use-histogram-matching" "0" "--winsorize-image-intensities" "[0.005,0.995]"
           "--initial-moving-transform" (str "[" target-mri "," source-mri ",1]")
           "--transform" "Affine[0.1]"
           "--metric" (str "MI[" target-mri "," source-mri ",1,32,Regular,0.25]")
           "--convergence" "[1000x500x250x100,1e-6,10]"
           "--shrink-factors" "8x4x2x1" "--smoothing-sigmas" "3x2x1x0vox")))

;;; HDbet
(defn hdbet
  "mask is saved as out_bet.nii.gz"
  [in out mask?]
  (shell `("hd-bet "
           "-i" ~in
           "-o" ~out
           ~(when mask? "--save_bet_mask")
           "--verbose")
         {:pre-start-fn (fn [v]
                          (println (clojure.string/join " " (:cmd v))))}))

;;; Core

(defn register
  [imgs]
  (let [;; Start by doing hdbet on all images (we keep the T1 mask)
        T1w-mask nil
        ;; Next register all masked images to the masked T1 image (keeping the registrations)
        
        ;; Next apply the registration to the unmasked images
        
        ;; Then apply the hdbet mask of the T1 image to the unmasked registered images
        
        ]))








;;;; LOGIC

;; (defn main []
;;   (let [source      (:s args)
;;         target      (:t args)
;;         participant (:p args)
;;         output-dir  (or (:d args) (System/getProperty "user.dir"))
;;         warp2mni?   (:w args)
;;         reg-type    (or (:r args) "1")]

;;     (when (and (nil? source) (nil? participant))
;;       (usage))

;;     (ensure-dir output-dir)

;;     ;; Direct registration
;;     (if source
;;       (let [source-label (-> (io/file source) .getName (str/split #"\.") first)
;;             target-label (-> (io/file target) .getName (str/split #"\.") first)
;;             out-prefix   (str output-dir "/" source-label "_reg2_" target-label)]

;;         (case reg-type
;;           "1" (rigid-register source target out-prefix)
;;           "2" (println "Affine registration not implemented yet")
;;           "3" (println "Non-rigid registration not implemented yet")))

;;       ;; Full participant mode (BIDS directory)
;;       (println "Participant mode is not fully implemented in BB yet."))))


