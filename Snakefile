rule all:
    input: "done.txt"

rule test:
    input: "one.txt"
    output: "done.txt"
    shell: "touch {output}"
