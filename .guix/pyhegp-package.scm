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

(define-module (pyhegp-package)
  #:use-module ((gnu packages check) #:select (python-hypothesis python-pytest))
  #:use-module ((gnu packages python-build) #:select (python-flit-core))
  #:use-module ((gnu packages python-science) #:select (python-pandas python-scipy))
  #:use-module ((gnu packages python-xyz) #:select (python-click python-numpy))
  #:use-module (guix build-system pyproject)
  #:use-module (guix gexp)
  #:use-module (guix git-download)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix packages)
  #:use-module (guix utils))

(define-public pyhegp
  (package
    (name "pyhegp")
    (version "0.1.0")
    (source (local-file ".."
                        "pyhegp-checkout"
                        #:recursive? #t
                        #:select? (lambda (file stat)
                                    ;; If .guix is included, changes
                                    ;; to other files under .guix—such
                                    ;; as the hsmice
                                    ;; test—unnecessarily trigger a
                                    ;; rebuild of pyhegp. This could
                                    ;; be a nuisance when hacking on
                                    ;; the test scripts.
                                    (and (not (string-contains file "/.guix/"))
                                         (not (string-contains file "/e2e-tests/"))
                                         ((or (git-predicate (dirname (current-source-directory)))
                                              (const #t))
                                          file stat)))))
    (build-system pyproject-build-system)
    (arguments
     (list #:phases
           #~(modify-phases %standard-phases
               ;; When the CI environment variable is set to any
               ;; value, hypothesis switches to using its built-in CI
               ;; profile.
               (add-before 'check 'set-ci-environment
                 (lambda _
                   (setenv "CI" "running-on-ci"))))))
    (inputs
     (list python-click
           python-numpy
           python-pandas
           python-scipy))
    (native-inputs
     (list python-flit-core
           python-hypothesis
           python-pytest))
    (home-page "https://github.com/encryption4genetics/pyhegp")
    (synopsis "Homomorphic encryption of genotypes and phenotypes")
    (description "@code{python-pyhegp} provides CLI utilities
implementing homomorphic encryption of genotypes and phenotypes.")
    (license license:gpl3+)))

pyhegp
