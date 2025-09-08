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

(define-module (readme-images)
  #:use-module ((gnu packages uml) #:select (plantuml))
  #:use-module (guix gexp))

(define readme-images-gexp
  (with-imported-modules '((guix build utils))
    #~(begin
        (use-modules (guix build utils))

        (invoke #$(file-append plantuml "/bin/plantuml")
                #$(local-file "../doc/simple-workflow.uml")
                #$(local-file "../doc/joint-workflow.uml")
                "-o" #$output))))

(define-public readme-images
  (computed-file "pyhegp-readme-images" readme-images-gexp))

readme-images
