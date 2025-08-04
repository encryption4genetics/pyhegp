#! /bin/sh

cat simple-workflow.uml | guix shell plantuml -- plantuml -p > simple-workflow.png
cat joint-workflow.uml | guix shell plantuml -- plantuml -p > joint-workflow.png
