(ns KUL-anat-register
  (:require [babashka.cli :as cli]
            [com.rpl.specter :as s]
            [util :refer [ensure-dir kulderivativesdir file-label
                          images-for-participant with-temp-file]]
            [register :refer [apply-transform
                              bids-register]]
            [tumor :refer [process-tumor]]
            [bet :refer [hdbet apply-mask]]
            [sc.api :refer [defsc spy]]))

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

(def interpolation-map {"1" "BSpline" "2" "NearestNeighbor" "3" "Linear"})
(def interpolation-type
  (get interpolation-map (str (:i args "1")) "BSpline"))
(def participant (:p args))
(def kul-log-dir (str "KUL_LOG/KUL_anat_register/sub-" participant))
(ensure-dir kul-log-dir)
;; (assert (get-bids-dir "."))
;; (def kulderivativesdir (str (get-bids-dir ".") "derivatives/KUL_compute/sub-" participant "/KUL_anat_register"))
;; (def registeroutputdir (str kulderivativesdir "/antsregister"))
(def verbose-level (or (:v args) 0))


(defn main
  [p dir]
  (let [imgs      (->> (images-for-participant p dir)
                       (s/transform [s/MAP-VALS]
                                    first))
        ;; Start by doing hdbet on all images (we keep the T1 mask)
        out-dir   (str (kulderivativesdir p dir) "/KUL_anat_register")
        _         (ensure-dir out-dir)
        tumor-dir (str (kulderivativesdir p dir) "/tumor")
        _         (ensure-dir tumor-dir)

        
        T1-image   (or (get imgs :T1w nil)
                       (throw (Exception. "No T1w image found")))
        other-imgs (dissoc imgs :T1w :SWI :SWIp :FGATIR)
        
        T1-bet  (str out-dir "/T1w_masked.nii.gz")
        T1-mask (str out-dir "/T1w_masked_bet.nii.gz")
        
        _ (hdbet T1-image T1-bet true)

        temp-hdbet     (fn [img]
                         (with-temp-file [temp-file ".nii.gz"]
                           (-> (hdbet img temp-file false) :out println)
                           {:native   img
                            :label    (file-label img)
                            :bet-file temp-file}))
        other-imgs-bet (s/transform [s/MAP-VALS] temp-hdbet other-imgs)
        
        ;; Next register all masked images to the masked T1 image (keeping the registrations)
        other-imgs-reg (s/transform
                        [s/MAP-VALS]
                        #(merge %
                                (bids-register {:output-dir   (str (kulderivativesdir p dir)
                                                                   "/KUL_anat_register")
                                                :source-mri   (:bet-file %)
                                                :type         :affine
                                                :label        (:label %)
                                                :target-mri   T1-bet
                                                :p            p
                                                :dir          dir
                                                :ants-verbose 0}))
                        other-imgs-bet)]
    
    ;; Next apply the registration to the unmasked images
    (run! #(let [in        (:native %)
                 output    (:output-mri %)
                 ref       (:T1w imgs)
                 transform (:warp-field %)]
             (apply-transform in output ref transform interpolation-type))
          (s/select [s/MAP-VALS] other-imgs-reg))
    
    ;; Then apply the hdbet mask of the T1 image to the unmasked registered images
    (run!
     #(apply-mask % T1-mask %)
     (s/select [s/MAP-VALS :output-mri]
               other-imgs-reg))
    ;; Copy all files to the processing dir (using a transformation to MNI space)
    (process-tumor p dir)
    ;; (let [{:keys [cT1w FLAIR T2w T1w]} (s/transform
    ;;                                     [s/MAP-VALS]
    ;;                                     :output-mri
    ;;                                     (assoc other-imgs-reg :T1w {:output-mri T1-bet}))
    
    ;;       mni-image  "/opt/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"
    ;;       warp-field (fs/create-temp-file)]
    
    ;;   (affine-register-command 0
    ;;                            warp-field
    ;;                            (str tumor-dir "/T1.nii.gz")
    ;;                            "BSpline"
    ;;                            mni-image
    ;;                            T1w)
    ;;   (let [transform (str warp-field "0GenericAffine.mat")]
    ;;     (apply-transform T2w (str tumor-dir "/T2.nii.gz") mni-image transform interpolation-type)
    ;;     (apply-transform cT1w (str tumor-dir "/CT1.nii.gz") mni-image transform interpolation-type)
    ;;     (apply-transform FLAIR (str tumor-dir "/FLAIR.nii.gz") mni-image transform interpolation-type))
    ;;   (run-command "hd_glio_predict" "-t1" (str tumor-dir "/T1.nii.gz")
    ;;                "-t1c" (str tumor-dir "/CT1.nii.gz")
    ;;                "-t2" (str tumor-dir "/T2.nii.gz")
    ;;                "-flair" (str tumor-dir "/FLAIR.nii.gz")
    ;;                "-o" (str tumor-dir "/tumor.nii.gz")))
    ))


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

