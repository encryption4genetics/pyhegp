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
  #:use-module ((gnu packages base) #:select (tar))
  #:use-module ((gnu packages compression) #:select (gzip))
  #:use-module ((gnu packages cran) #:select (r-dplyr r-purrr r-qqman r-stringr))
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

(define test-profile
  (profile
   (content (packages->manifest (list gzip tar pyhegp
                                      python python-click python-pandas
                                      r r-dplyr r-genio
                                      r-mixed-model-gwas r-purrr
                                      r-qqman r-readr r-stringr
                                      r-tibble r-tidyr)))))

(define wrangle-script
  (local-file "../e2e-tests/hsmice/wrangle.r"))

(define gwas-script
  (local-file "../e2e-tests/hsmice/gwas.r"))

(define check-qtl-script
  (local-file "../e2e-tests/hsmice/check-qtl.py"))

(define hsmice-test-gexp
  (with-imported-modules '((guix build utils))
    #~(begin
        (use-modules (guix build utils))

        (mkdir #$output)
        (set-path-environment-variable
         "PATH" '("/bin") '(#$test-profile))
        (set-path-environment-variable
         "GUIX_PYTHONPATH"
         '(#$(string-append "/lib/python"
                            (version-major+minor (package-version python))
                            "/site-packages"))
         '(#$test-profile))
        (set-path-environment-variable
         "R_LIBS_SITE" '("/site-library") '(#$test-profile))
        (invoke "tar" "-xvf" #$hsmice-data
                "./HSmice/1_QTL_data/")
        (invoke "Rscript" #$wrangle-script "HSmice/1_QTL_data" ".")

        ;; GWAS on plaintext
        (invoke "Rscript" #$gwas-script
                "genotype.tsv" "phenotype.tsv"
                (string-append #$output "/plaintext-pvalues"))
        (copy-file "Rplots.pdf" (string-append #$output "/plaintext-manhattan.pdf"))

        ;; GWAS with simple ciphertext data sharing
        (invoke "pyhegp" "encrypt" "genotype.tsv" "phenotype.tsv")
        (invoke "Rscript" #$gwas-script
                "genotype.tsv.hegp" "phenotype.tsv.hegp"
                (string-append #$output "/ciphertext-pvalues"))
        (copy-file "Rplots.pdf"
                   (string-append #$output "/ciphertext-manhattan.pdf"))

        ;; Joint federated GWAS
        (invoke "pyhegp" "summary" "genotype1.tsv" "-o" "summary1")
        (invoke "pyhegp" "summary" "genotype2.tsv" "-o" "summary2")
        (invoke "pyhegp" "pool" "-o" "complete-summary" "summary1" "summary2")
        (invoke "pyhegp" "encrypt" "-s" "complete-summary" "genotype1.tsv" "phenotype1.tsv")
        (invoke "pyhegp" "encrypt" "-s" "complete-summary" "genotype2.tsv" "phenotype2.tsv")
        (invoke "pyhegp" "cat-genotype" "-o" "complete-genotype.tsv.hegp"
                "genotype1.tsv.hegp" "genotype2.tsv.hegp")
        (invoke "pyhegp" "cat-phenotype" "-o" "complete-phenotype.tsv.hegp"
                "phenotype1.tsv.hegp" "phenotype2.tsv.hegp")
        (invoke "Rscript" #$gwas-script
                "complete-genotype.tsv.hegp" "complete-phenotype.tsv.hegp"
                (string-append #$output "/federated-ciphertext-pvalues"))
        (copy-file "Rplots.pdf"
                   (string-append #$output "/federated-ciphertext-manhattan.pdf"))

        ;; Check that the QTL is where it should be.
        (for-each (lambda (pvalues-file)
                    (invoke "python3" #$check-qtl-script
                            (string-append #$output "/" pvalues-file)))
                  (list "plaintext-pvalues"
                        "ciphertext-pvalues"
                        "federated-ciphertext-pvalues")))))

(define-public hsmice-test
  (computed-file "hsmice-test" hsmice-test-gexp))

hsmice-test
