import yaml, sys

cfg = yaml.safe_load(open(sys.argv[1], "r", encoding="utf-8"))
out = sys.argv[2]

with open(out, "w", encoding="utf-8") as f:
    f.write("#pragma once\n")
    f.write("namespace CommonParams {\n")
    f.write(f"  inline constexpr int adc_peak_min_index = {int(cfg['adc_peak_min_index'])};\n")
    f.write(f"  inline constexpr int adc_peak_max_index = {int(cfg['adc_peak_max_index'])};\n")
    f.write(f"  inline constexpr int example_channel = {int(cfg['example_channel'])};\n")
    f.write("}\n")