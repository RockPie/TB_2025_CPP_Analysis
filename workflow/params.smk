rule gen_params_header:
    input:
        "config/params.yaml",
        "workflow/gen_params_header.py"
    output:
        "include/CommonParams.hxx"
    shell:
        "python3 {input[1]} {input[0]} {output}"