#! /bin/sh

cat workflow.uml | guix shell plantuml -- plantuml -p > workflow.png
