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

(use-modules ((gnu packages check) #:select (python-pytest))
             ((gnu packages python-build) #:select (python-flit-core))
             ((gnu packages python-science) #:select (python-scipy))
             ((gnu packages python-xyz) #:select (python-click python-numpy))
             (guix build-system pyproject)
             (guix gexp)
             (guix git-download)
             ((guix licenses) #:prefix license:)
             (guix packages)
             (guix utils))

(define-public python-pyhegp
  (package
    (name "python-pyhegp")
    (version "0.1.0")
    (source (local-file "."
                        "pyhegp-checkout"
                        #:recursive? #t
                        #:select? (git-predicate (current-source-directory))))
    (build-system pyproject-build-system)
    (arguments
     (list #:tests? #f))
    (native-inputs
     (list python-flit-core
           python-pytest))
    (propagated-inputs
     (list python-click
           python-numpy
           python-scipy))
    (home-page "https://github.com/encryption4genetics/pyhegp")
    (synopsis "Homomorphic encryption of genotypes and phenotypes")
    (description "@code{python-pyhegp} provides a Python library and CLI utilities
implementing homomorphic encryption of genotypes and phenotypes.")
    (license license:gpl3+)))

python-pyhegp
