;;; pyhegp --- Homomorphic encryption of genotypes and phenotypes
;;; Copyright © 2025 Arun Isaac <arunisaac@systemreboot.net>
;;;
;;; This file is part of pyhegp.
;;;
;;; pyhegp is free software: you can redistribute it and/or modify it
;;; under the terms of the GNU General Public License as published by
;;; the Free Software Foundation, either version 3 of the License, or
;;; (at your option) any later version.
;;;
;;; pyhegp is distributed in the hope that it will be useful, but
;;; WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
;;; General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with pyhegp. If not, see <https://www.gnu.org/licenses/>.

(define-module (hsmice-test)
  #:use-module ((gn packages bioinformatics) #:select (r-genio))
  #:use-module ((gn packages julia) #:select (julia-pipe))
  #:use-module ((gn packages julia) #:select (julia-jwas) #:prefix guix:)
  #:use-module ((gnu packages base) #:select (tar))
  #:use-module ((gnu packages compression) #:select (gzip))
  #:use-module ((gnu packages cran) #:select (r-dplyr r-purrr r-qqman r-stringr))
  #:use-module ((gnu packages julia) #:select (julia))
  #:use-module ((gnu packages julia-xyz) #:select (julia-csv julia-dataframes))
  #:use-module ((gnu packages python) #:select (python))
  #:use-module ((gnu packages python-science) #:select (python-pandas))
  #:use-module ((gnu packages python-xyz) #:select (python-click))
  #:use-module ((gnu packages statistics) #:select (r r-readr r-tibble r-tidyr))
  #:use-module (guix build-system r)
  #:use-module (guix download)
  #:use-module (guix gexp)
  #:use-module (guix git-download)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix packages)
  #:use-module (guix profiles)
  #:use-module (guix utils)
  #:use-module ((pyhegp-package) #:select (pyhegp)))

;; See details of dataset at
;; https://rdr.ucl.ac.uk/articles/dataset/HSmice_tar_gz/24114471?file=42304248
(define hsmice-data
  (origin
    (method url-fetch)
    (uri "https://ndownloader.figshare.com/files/42304248")
    (file-name "HSmice.tar.gz")
    (sha256
     (base32
      "1s6a83r0mll8z2lfv1b94zr2sjdrky5nyq1mpgl8fjjb5s8v2vyx"))))

(define-public r-mixed-model-gwas
  (package
   (name "r-mixed-model-gwas")
   (version "1.3.1")
   (source (origin
            (method git-fetch)
            (uri (git-reference
                  (url "https://github.com/encryption4genetics/mixed-model-gwas")
                  (commit (string-append "v" version))))
            (file-name (git-file-name name version))
            (sha256
             (base32
              "0vll55v8wjc0179n5q9ch9ah3dvgymc374wlbz33yzyi35yr8ds2"))))
   (build-system r-build-system)
   (home-page "https://github.com/encryption4genetics/mixed-model-gwas")
   (synopsis "R mixed model GWAS")
   (description "@code{r-mixed-model-gwas} implements a mixed model @acronym{GWAS,
genome-wide association study} library for R.")
   (license license:gpl3+)))

(define hsmice-wrangled-gexp
  (let ((script-profile (profile
                         (content (packages->manifest
                                   (list gzip tar r r-dplyr r-genio
                                         r-purrr r-readr r-tibble r-tidyr))))))
    (with-imported-modules '((guix build utils))
      #~(begin
          (use-modules (guix build utils))

          (mkdir #$output)
          (set-path-environment-variable
           "PATH" '("bin") '(#$script-profile))
          (set-path-environment-variable
           "R_LIBS_SITE" '("site-library") '(#$script-profile))
          (invoke "tar" "-xvf" #$hsmice-data
                  "./HSmice/1_QTL_data/")
          (invoke "Rscript"
                  #$(local-file "../e2e-tests/hsmice/wrangle.r")
                  "HSmice/1_QTL_data" #$output)))))

(define hsmice-wrangled
  (computed-file "hsmice-wrangled" hsmice-wrangled-gexp))

(define hsmice-ciphertext-gexp
  (let ((script-profile (profile
                          (content (packages->manifest (list pyhegp))))))
    (with-imported-modules '((guix build utils))
      #~(begin
          (use-modules (guix build utils)
                       (srfi srfi-26))

          (mkdir #$output)
          (set-path-environment-variable
           "PATH" '("bin") '(#$script-profile))
          (for-each (cut install-file <> (getcwd))
                    (find-files #$hsmice-wrangled "\\.tsv$"))
          ;; Simple data sharing workflow
          (invoke "pyhegp" "encrypt" "genotype.tsv" "phenotype.tsv")
          ;; Joint/federated analysis workflow
          (invoke "pyhegp" "summary" "genotype1.tsv" "-o" "summary1")
          (invoke "pyhegp" "summary" "genotype2.tsv" "-o" "summary2")
          (invoke "pyhegp" "pool" "-o" "complete-summary" "summary1" "summary2")
          (invoke "pyhegp" "encrypt" "-s" "complete-summary" "genotype1.tsv" "phenotype1.tsv")
          (invoke "pyhegp" "encrypt" "-s" "complete-summary" "genotype2.tsv" "phenotype2.tsv")
          (invoke "pyhegp" "cat-genotype" "-o" "complete-genotype.tsv.hegp"
                  "genotype1.tsv.hegp" "genotype2.tsv.hegp")
          (invoke "pyhegp" "cat-phenotype" "-o" "complete-phenotype.tsv.hegp"
                  "phenotype1.tsv.hegp" "phenotype2.tsv.hegp")
          (for-each (cut install-file <> #$output)
                    (find-files (getcwd) "\\.tsv.hegp$"))))))

(define hsmice-ciphertext
  (computed-file "hsmice-ciphertext" hsmice-ciphertext-gexp))

(define hsmice-r-mixed-model-gwas-gexp
  (let ((gwas-script (local-file "../e2e-tests/hsmice/gwas.r"))
        (script-profile (profile
                          (content (packages->manifest
                                    (list r r-dplyr r-mixed-model-gwas
                                          r-qqman r-readr r-stringr
                                          r-tibble r-tidyr))))))
    (with-imported-modules '((guix build utils))
      #~(begin
          (use-modules (guix build utils))

          (mkdir #$output)
          (set-path-environment-variable
           "PATH" '("bin") '(#$script-profile))
          (set-path-environment-variable
           "R_LIBS_SITE" '("site-library") '(#$script-profile))

          ;; GWAS on plaintext
          (invoke "Rscript" #$gwas-script
                  #$(file-append hsmice-wrangled "/genotype.tsv")
                  #$(file-append hsmice-wrangled "/phenotype.tsv")
                  (string-append #$output "/plaintext-pvalues"))
          (copy-file "Rplots.pdf"
                     (string-append #$output "/plaintext-manhattan.pdf"))

          ;; GWAS with simple ciphertext data sharing
          (invoke "Rscript" #$gwas-script
                  #$(file-append hsmice-ciphertext "/genotype.tsv.hegp")
                  #$(file-append hsmice-ciphertext "/phenotype.tsv.hegp")
                  (string-append #$output "/ciphertext-pvalues"))
          (copy-file "Rplots.pdf"
                     (string-append #$output "/ciphertext-manhattan.pdf"))

          ;; Joint federated GWAS
          (invoke "Rscript" #$gwas-script
                  #$(file-append hsmice-ciphertext "/complete-genotype.tsv.hegp")
                  #$(file-append hsmice-ciphertext "/complete-phenotype.tsv.hegp")
                  (string-append #$output "/federated-ciphertext-pvalues"))
          (copy-file "Rplots.pdf"
                     (string-append #$output "/federated-ciphertext-manhattan.pdf"))))))

(define hsmice-r-mixed-model-gwas
  (computed-file "hsmice-r-mixed-model-gwas" hsmice-r-mixed-model-gwas-gexp))

(define hsmice-qtl-checked-gexp
  (let ((script-profile (profile
                          (content (packages->manifest
                                    (list python python-pandas))))))
    (with-imported-modules '((guix build utils))
      #~(begin
          (use-modules (guix build utils)
                       (srfi srfi-26))

          (mkdir #$output)
          (set-path-environment-variable
           "PATH" '("bin") '(#$script-profile))
          (set-path-environment-variable
           "GUIX_PYTHONPATH"
           '(#$(string-append "lib/python"
                              (version-major+minor (package-version python))
                              "/site-packages"))
           '(#$script-profile))

          ;; Check that the QTL is where it should be.
          (for-each (cut invoke
                         "python3"
                         #$(local-file "../e2e-tests/hsmice/check-qtl.py")
                         <>
			 "p < 1e-10")
                    (find-files #$hsmice-r-mixed-model-gwas
                                "\\-pvalues$"))))))

(define hsmice-qtl-checked
  (computed-file "hsmice-qtl-checked" hsmice-qtl-checked-gexp))

;; Our CI system does not have enough memory to run tests for
;; julia-jwas. So, we disable them.
(define julia-jwas
  (package
    (inherit guix:julia-jwas)
    (arguments
     (substitute-keyword-arguments (package-arguments guix:julia-jwas)
       ((#:tests? tests? #f) #f)))))

(define hsmice-jwas-gwas-gexp
  (let ((gwas-script (local-file "../e2e-tests/hsmice/jwas-gwas.jl"))
        (jwas-manhattan-script (local-file "../e2e-tests/hsmice/jwas-manhattan.r"))
        (script-profile (profile
                          (content (packages->manifest
				    (list julia julia-csv julia-dataframes
					  julia-jwas julia-pipe
                                          r r-dplyr r-qqman r-readr))))))
    (with-imported-modules '((guix build utils))
      #~(begin
          (use-modules (guix build utils))

          (mkdir #$output)
	  (setenv "HOME" "/tmp")
          (set-path-environment-variable
           "PATH" '("bin") '(#$script-profile))
          (set-path-environment-variable
           "JULIA_LOAD_PATH" '("share/julia/loadpath") '(#$script-profile))
	  (set-path-environment-variable
           "JULIA_DEPOT_PATH" '("share/julia") '(#$script-profile))
          (set-path-environment-variable
           "R_LIBS_SITE" '("site-library") '(#$script-profile))

          ;; GWAS on plaintext
          (invoke "julia" #$gwas-script
                  #$(file-append hsmice-wrangled "/genotype.tsv")
                  #$(file-append hsmice-wrangled "/phenotype.tsv")
                  (string-append #$output "/plaintext-jwas-gwas.tsv"))
          (invoke "Rscript" #$jwas-manhattan-script
                  (string-append #$output "/plaintext-jwas-gwas.tsv"))
          (copy-file "Rplots.pdf"
                     (string-append #$output "/plaintext-jwas-manhattan.pdf"))

          ;; GWAS with simple ciphertext data sharing
          (invoke "julia" #$gwas-script
                  #$(file-append hsmice-ciphertext "/genotype.tsv.hegp")
                  #$(file-append hsmice-ciphertext "/phenotype.tsv.hegp")
                  (string-append #$output "/ciphertext-jwas-gwas.tsv"))
          (invoke "Rscript" #$jwas-manhattan-script
                  (string-append #$output "/ciphertext-jwas-gwas.tsv"))
          (copy-file "Rplots.pdf"
                     (string-append #$output "/ciphertext-jwas-manhattan.pdf"))

          ;; Joint federated GWAS
          (invoke "julia" #$gwas-script
                  #$(file-append hsmice-ciphertext "/complete-genotype.tsv.hegp")
                  #$(file-append hsmice-ciphertext "/complete-phenotype.tsv.hegp")
                  (string-append #$output "/federated-ciphertext-jwas-gwas.tsv"))
          (invoke "Rscript" #$jwas-manhattan-script
                  (string-append #$output "/federated-ciphertext-jwas-gwas.tsv"))
          (copy-file "Rplots.pdf"
                     (string-append #$output "/federated-ciphertext-jwas-manhattan.pdf"))))))

(define hsmice-jwas-gwas
  (computed-file "hsmice-jwas-gwas" hsmice-jwas-gwas-gexp))

(define hsmice-jwas-qtl-checked-gexp
  (let ((script-profile (profile
                          (content (packages->manifest
                                    (list python python-pandas))))))
    (with-imported-modules '((guix build utils))
      #~(begin
          (use-modules (guix build utils)
                       (srfi srfi-26))

          (mkdir #$output)
          (set-path-environment-variable
           "PATH" '("bin") '(#$script-profile))
          (set-path-environment-variable
           "GUIX_PYTHONPATH"
           '(#$(string-append "lib/python"
                              (version-major+minor (package-version python))
                              "/site-packages"))
           '(#$script-profile))

          ;; Check that the QTL is where it should be.
          (for-each (cut invoke
                         "python3"
                         #$(local-file "../e2e-tests/hsmice/check-qtl.py")
                         <>
			 "modelfrequency > 0.4")
                    (find-files #$hsmice-jwas-gwas
                                "\\-gwas.tsv$"))))))

(define hsmice-jwas-qtl-checked
  (computed-file "hsmice-jwas-qtl-checked" hsmice-jwas-qtl-checked-gexp))

(define-public hsmice-test
  (directory-union "hsmice-test"
                   (list hsmice-r-mixed-model-gwas
                         hsmice-qtl-checked
			 hsmice-jwas-gwas
			 hsmice-jwas-qtl-checked)))

hsmice-test
