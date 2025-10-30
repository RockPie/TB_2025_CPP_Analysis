import lib_plotting
import ROOT, os, json, time
import numpy as np

ROOT.gROOT.SetBatch(True)

# load color configuration
color_json_file = "config/run_color_collection.json"
color_json = json.load(open(color_json_file, "r"))
channel_color_map = {}
for color_hex, run_info in color_json.items():
    run_list = run_info.get("run_numbers", [])
    r = int(color_hex[1:3], 16)
    g = int(color_hex[3:5], 16)
    b = int(color_hex[5:7], 16)
    color_code = ROOT.TColor.GetColor(r, g, b)
    for run_number in run_list:
        channel_color_map[run_number] = {}
        channel_color_map[run_number]["color"] = color_code
    # add description printout
    description = run_info.get("description", "")
    channel_color_map[run_number]["description"] = description

list_root_file_folders = [
    "dump/301_EventReconX/beamtests",
    "dump/302_EventMatchX/beamtests",
    "dump/303_Pedestal/beamtests"
]
list_root_file_info_str = [
    "EventReconX_Beamtests",
    "EventMatchX_Beamtests",
    "Pedestal_Beamtests"
]



# for each folder, create plots of all TParameters
for folder_index, root_file_folder in enumerate(list_root_file_folders):
    set_TParameter_names = set()
    list_TParameter_lists = []
    list_TParameter_runs = []
    root_file_list = [f for f in os.listdir(root_file_folder) if f.endswith(".root")]
    run_numbers = []
    for root_file_name in root_file_list:
        parts = root_file_name.split("_")
        for part in parts:
            if part.startswith("Run"):
                run_number_str = part.replace("Run", "")
                run_number_str = run_number_str.replace(".root", "")
                if run_number_str.isdigit():
                    run_numbers.append(int(run_number_str))
    print(f"Found {len(run_numbers)} runs in folder {root_file_folder}")
    print(f"Run numbers: {run_numbers}")

    # go through each root file to figure out TParameters
    for root_file_name in root_file_list:
        root_file_path = os.path.join(root_file_folder, root_file_name)
        # print(f"Processing file: {root_file_path}")
        root_file = ROOT.TFile.Open(root_file_path, "READ")
        if not root_file or root_file.IsZombie():
            print(f"Error opening file: {root_file_path}")
            continue

        list_keys = root_file.GetListOfKeys()
        for key in list_keys:
            key_name = key.GetName()
            key_class_name = key.GetClassName()
            if key_class_name != "TParameter<double>" and key_class_name != "TParameter<int>":
                continue
            if key_name not in set_TParameter_names:
                print(f"Found new TParameter: {key_name} of type {key_class_name}")
                index = -1
                TParameter_obj = root_file.Get(key_name)
                if not TParameter_obj:
                    print(f"Error retrieving TParameter {key_name} from file {root_file_path}")
                    continue
                set_TParameter_names.add(key_name)
                list_TParameter_lists.append([])
                list_TParameter_lists[-1].append(TParameter_obj.GetVal())
                list_TParameter_runs.append(run_number_str)
            else:
                index = list(set_TParameter_names).index(key_name)
                TParameter_obj = root_file.Get(key_name)
                if TParameter_obj:
                    TParameter_value = TParameter_obj.GetVal()
                    list_TParameter_lists[index].append(TParameter_value)
                    list_TParameter_runs.append(run_number_str)
        root_file.Close()

    # create plots for each TParameter
    output_plot_folder = "qa_plots/TParameter_PlotsX/" + list_root_file_info_str[folder_index]
    os.makedirs(output_plot_folder, exist_ok=True)

    # sort the values in each TParameter list by run number
    for index in range(len(list_TParameter_lists)):
        list_TParameter_lists[index].sort(key=lambda x: run_numbers.index(int(list_TParameter_runs[index])))
    sorted_run_numbers = sorted(run_numbers)
    run_to_color =  {r: info["color"] for r, info in channel_color_map.items()}
    run_to_text  = {r: info.get("description", "") for r, info in channel_color_map.items()}

    for index, TParameter_name in enumerate(set_TParameter_names):
        title_lines = [
            "FoCal-H Prototype 3",
            "Beam Test 2025 Oct",
            f"TParameter: {TParameter_name}",
            time.strftime("Date: %Y-%m-%d")
        ]
        out_pdf_file = os.path.join(output_plot_folder, f"TParameter_{TParameter_name}.pdf")

        # get the 90% max y
        y_max_90 = np.percentile(list_TParameter_lists[index], 90)

        lib_plotting.draw_run_categorical_bands_scatter(
            runs=sorted_run_numbers,
            y_values=[x for x in list_TParameter_lists[index]],
            y_errors=[0 for _ in list_TParameter_lists[index]],
            run_to_color=run_to_color,
            output_pdf=out_pdf_file,
            y_min=0,
            y_max=y_max_90 * 1.5,
            logy=False,
            band_shrink=0.00,
            y_label=TParameter_name,
            x_label="Run Number",
            title_lines=title_lines,
            run_to_text=run_to_text,
            png_also=False,
            root_also=os.path.join(output_plot_folder, f"TParameter_{TParameter_name}.root")
        )

