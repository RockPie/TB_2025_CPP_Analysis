import ROOT, os, json, time, math

def draw_run_categorical_bands_multi(
    runs,                    # list[int]            —— 原始 run 列表（可重复）
    y_value_lists,                # list[float]          —— 与 runs 等长
    y_error_lists,          # list[float] or None  —— 与 runs 等长，若提供则画误差棒
    run_to_color,            # dict[int]->int       —— run号 -> ROOT 颜色索引
    output_pdf,              # str                  —— 输出 PDF 路径
    y_min=1e-6,              # float                —— y 轴最小值（logY 必须 > 0）
    y_max=1.0,               # float
    logy=True,               # bool                 —— 是否 logY
    band_shrink=0.02,        # float                —— 每个色带左右各收缩的比例
    y_label="Y",             # str                  —— y 轴标题
    x_label="Run Number",    # str                  —— x 轴标题
    title_lines=None,        # list[str]            —— 左上多行标题
    run_to_text=None,        # dict[int]->str       —— “仅第一次出现时”在色带内竖直标注（可选）
    png_also=True,           # bool                 —— 额外导出 PNG
    root_also=None,          # Optional[str]        —— 若提供则把 canvas 写入 ROOT 文件
    rotate_x_labels=True,    # bool                 —— x 轴标签竖排
    marker_style=20,         # int
    marker_size=0.6,         # float
    marker_colors=None # list[int]
):
    y_values = y_value_lists[0]
    y_errors = y_error_lists[0] if y_error_lists is not None else None
    assert len(runs) == len(y_values), "runs mismatch y_values length"
    if y_errors is not None:
        assert len(runs) == len(y_errors), "runs mismatch y_errors length"

    # 过滤/夹值（logY 不允许 0 或负数）
    vals = [y_min if (v is None or (isinstance(v, float) and math.isnan(v)) or v <= 0) else float(v)
            for v in y_values]
    if y_errors is not None:
        for i in range(len(vals)):
            err = y_errors[i]
            if err < 0:
                err = abs(err)
            y_lo = vals[i] - err
            if logy and y_lo <= 0:
                err = vals[i] - y_min
            vals[i] = vals[i]
            y_errors[i] = err
    runs = list(map(int, runs))

    # 分类轴顺序：按出现顺序去重
    runs_seen, seen = [], set()
    for r in runs:
        if r not in seen:
            runs_seen.append(r); seen.add(r)
    N = len(runs_seen)
    if N == 0:
        raise RuntimeError("空的 run 列表")

    # 画布
    c = ROOT.TCanvas(f"c_{int(time.time()*1000)}", "QA", 1600, 560)
    c.SetLeftMargin(0.10)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.18)
    c.SetTopMargin(0.10)
    c.SetLogy(bool(logy))

    # 框架：N 个等宽 bin，x∈[0, N]
    frame = ROOT.TH2F(f"frame_{int(time.time()*1000)}",
                      f";{x_label};{y_label}",
                      N, 0.0, float(N), 10, float(y_min), float(y_max))
    frame.SetStats(0); frame.SetFillStyle(0); frame.SetLineColor(0)
    for i, r in enumerate(runs_seen, start=1):
        frame.GetXaxis().SetBinLabel(i, str(r))
    if rotate_x_labels:
        frame.GetXaxis().LabelsOption("v")
    frame.GetXaxis().SetLabelSize(0.030)
    frame.GetYaxis().SetTitleOffset(0.7)
    frame.GetYaxis().SetLabelSize(0.030)
    frame.Draw("AXIS")

    # 背景色带（保存引用避免 GC）——仅在下一列不一样时写
    boxes, texts, titles = [], [], []
    for i, r in enumerate(runs_seen):
        x1, x2 = float(i), float(i + 1)
        if band_shrink > 0:
            s = band_shrink * (x2 - x1)
            x1 += s; x2 -= s

        col = run_to_color.get(r, ROOT.kGray+1)
        box = ROOT.TBox(x1, float(y_min), x2, float(y_max))
        box.SetLineColor(0)
        box.SetFillStyle(1001)   # 实心，不透明 —— PDF 最稳
        box.SetFillColor(int(col))
        box.Draw("same f")
        boxes.append(box)

        if run_to_text:
            desc = run_to_text.get(r, "")
            r_next = runs_seen[i + 1] if i < len(runs_seen) - 1 else None
            desc_next = run_to_text.get(r_next, "") if r_next is not None else ""
            if desc_next == "" or desc != desc_next:
                # skip if this is the left most column to avoid overlap with y axis
                # if i <= 1:
                #     continue
                # printed_desc.add(desc)
                t = ROOT.TLatex()
                t.SetTextFont(52)
                t.SetTextAlign(13)
                t.SetTextSize(0.025)
                t.SetTextColor(ROOT.kGray+2)
                t.SetTextAngle(90)
                if not logy:
                    y_0_9 = 0.05
                else: # calculate the 0.9 position in log scale
                    y_0_9 = math.exp(math.log(y_min) + 0.05 * (math.log(y_max) - math.log(y_min)))
                t.DrawLatex(x1 + 0.2*(x2 - x1), y_0_9, str(desc))
                texts.append(t)

    # 散点（画在 bin 中心）
    pos_map = {r: i for i, r in enumerate(runs_seen)}
    for idx, (y_values, y_errors) in enumerate(zip(y_value_lists, y_error_lists)):
        g = ROOT.TGraphErrors(len(runs)) if y_errors is not None else ROOT.TGraph(len(runs))
        if y_errors is not None:
            for i, err in enumerate(y_errors):
                r, y = runs[i], y_values[i]
                x = pos_map[r] + 0.5
                g.SetPoint(i, float(x), float(y))
                g.SetPointError(i, 0.0, float(err))
        else:
            for i, (r, y) in enumerate(zip(runs, y_values)):
                x = pos_map[r] + 0.5
                g.SetPoint(i, float(x), float(y))
        g.SetMarkerStyle(marker_style)
        g.SetMarkerSize(marker_size)
        if marker_colors and idx < len(marker_colors):
            g.SetMarkerColor(marker_colors[idx])
        else:
            g.SetMarkerColor(ROOT.kBlack)
        g.Draw("P SAME")
    # g = ROOT.TGraphErrors(len(runs)) if y_errors is not None else ROOT.TGraph(len(runs))
    # if y_errors is not None:
    #     for i, err in enumerate(y_errors):
    #         r, y = runs[i], vals[i]
    #         x = pos_map[r] + 0.5
    #         g.SetPoint(i, float(x), float(y))
    #         g.SetPointError(i, 0.0, float(err))
    # else:
    #     for i, (r, y) in enumerate(zip(runs, vals)):
    #         x = pos_map[r] + 0.5
    #         g.SetPoint(i, float(x), float(y))
    # g.SetMarkerStyle(marker_style)
    # g.SetMarkerSize(marker_size)
    # g.SetMarkerColor(marker_color)
    # g.Draw("P SAME")

    if title_lines:
        t0 = ROOT.TLatex()
        t0.SetNDC(True)
        y0 = 0.85
        for j, line in enumerate(title_lines):
            t0.SetTextFont(62 if j == 0 else 42)
            t0.SetTextSize(0.045 if j == 0 else 0.035)
            t0.DrawLatex(0.12, y0 - 0.045*j, line)
        titles.append(t0)

    ROOT.gPad.RedrawAxis()
    c.Update()

    # 输出
    os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
    c.SaveAs(output_pdf)
    if png_also:
        c.SaveAs(output_pdf.replace(".pdf", ".png"))
    if root_also:
        rfdir = os.path.dirname(root_also)
        if rfdir: os.makedirs(rfdir, exist_ok=True)
        rf = ROOT.TFile.Open(root_also, "RECREATE")
        c.Write("canvas")
        rf.Close()

    # 返回对象，避免被 GC
    return {
        "canvas": c,
        "frame": frame,
        "graph": g,
        "boxes": boxes,
        "texts": texts,
        "titles": titles,
        "runs_seen": runs_seen
    }
    

def draw_run_categorical_bands_scatter(
    runs,                    # list[int]            —— 原始 run 列表（可重复）
    y_values,                # list[float]          —— 与 runs 等长
    y_errors,          # list[float] or None  —— 与 runs 等长，若提供则画误差棒
    run_to_color,            # dict[int]->int       —— run号 -> ROOT 颜色索引
    output_pdf,              # str                  —— 输出 PDF 路径
    y_min=1e-6,              # float                —— y 轴最小值（logY 必须 > 0）
    y_max=1.0,               # float
    logy=True,               # bool                 —— 是否 logY
    band_shrink=0.02,        # float                —— 每个色带左右各收缩的比例
    y_label="Y",             # str                  —— y 轴标题
    x_label="Run Number",    # str                  —— x 轴标题
    title_lines=None,        # list[str]            —— 左上多行标题
    run_to_text=None,        # dict[int]->str       —— “仅第一次出现时”在色带内竖直标注（可选）
    png_also=True,           # bool                 —— 额外导出 PNG
    root_also=None,          # Optional[str]        —— 若提供则把 canvas 写入 ROOT 文件
    rotate_x_labels=True,    # bool                 —— x 轴标签竖排
    marker_style=20,         # int
    marker_size=0.6,         # float
    marker_color=ROOT.kBlack # int
):
    assert len(runs) == len(y_values), "runs mismatch y_values length"
    if y_errors is not None:
        assert len(runs) == len(y_errors), "runs mismatch y_errors length"

    # 过滤/夹值（logY 不允许 0 或负数）
    vals = [y_min if (v is None or (isinstance(v, float) and math.isnan(v)) or v <= 0) else float(v)
            for v in y_values]
    if y_errors is not None:
        for i in range(len(vals)):
            err = y_errors[i]
            if err < 0:
                err = abs(err)
            y_lo = vals[i] - err
            if logy and y_lo <= 0:
                err = vals[i] - y_min
            vals[i] = vals[i]
            y_errors[i] = err
    runs = list(map(int, runs))

    # 分类轴顺序：按出现顺序去重
    runs_seen, seen = [], set()
    for r in runs:
        if r not in seen:
            runs_seen.append(r); seen.add(r)
    N = len(runs_seen)
    if N == 0:
        raise RuntimeError("空的 run 列表")

    # 画布
    c = ROOT.TCanvas(f"c_{int(time.time()*1000)}", "QA", 1600, 560)
    c.SetLeftMargin(0.10)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.18)
    c.SetTopMargin(0.10)
    c.SetLogy(bool(logy))

    # 框架：N 个等宽 bin，x∈[0, N]
    frame = ROOT.TH2F(f"frame_{int(time.time()*1000)}",
                      f";{x_label};{y_label}",
                      N, 0.0, float(N), 10, float(y_min), float(y_max))
    frame.SetStats(0); frame.SetFillStyle(0); frame.SetLineColor(0)
    for i, r in enumerate(runs_seen, start=1):
        frame.GetXaxis().SetBinLabel(i, str(r))
    if rotate_x_labels:
        frame.GetXaxis().LabelsOption("v")
    frame.GetXaxis().SetLabelSize(0.030)
    frame.GetYaxis().SetTitleOffset(0.7)
    frame.GetYaxis().SetLabelSize(0.030)
    frame.Draw("AXIS")

    # 背景色带（保存引用避免 GC）——仅在下一列不一样时写
    boxes, texts, titles = [], [], []
    for i, r in enumerate(runs_seen):
        x1, x2 = float(i), float(i + 1)
        if band_shrink > 0:
            s = band_shrink * (x2 - x1)
            x1 += s; x2 -= s

        col = run_to_color.get(r, ROOT.kGray+1)
        box = ROOT.TBox(x1, float(y_min), x2, float(y_max))
        box.SetLineColor(0)
        box.SetFillStyle(1001)   # 实心，不透明 —— PDF 最稳
        box.SetFillColor(int(col))
        box.Draw("same f")
        boxes.append(box)

        if run_to_text:
            desc = run_to_text.get(r, "")
            r_next = runs_seen[i + 1] if i < len(runs_seen) - 1 else None
            desc_next = run_to_text.get(r_next, "") if r_next is not None else ""
            if desc_next == "" or desc != desc_next:
                # skip if this is the left most column to avoid overlap with y axis
                # if i <= 1:
                #     continue
                # printed_desc.add(desc)
                t = ROOT.TLatex()
                t.SetTextFont(52)
                t.SetTextAlign(13)
                t.SetTextSize(0.025)
                t.SetTextColor(ROOT.kGray+2)
                t.SetTextAngle(90)
                if not logy:
                    y_0_9 = 0.05
                else: # calculate the 0.9 position in log scale
                    y_0_9 = math.exp(math.log(y_min) + 0.05 * (math.log(y_max) - math.log(y_min)))
                t.DrawLatex(x1 + 0.2*(x2 - x1), y_0_9, str(desc))
                texts.append(t)

    # 散点（画在 bin 中心）
    pos_map = {r: i for i, r in enumerate(runs_seen)}
    g = ROOT.TGraphErrors(len(runs)) if y_errors is not None else ROOT.TGraph(len(runs))
    if y_errors is not None:
        for i, err in enumerate(y_errors):
            r, y = runs[i], vals[i]
            x = pos_map[r] + 0.5
            g.SetPoint(i, float(x), float(y))
            g.SetPointError(i, 0.0, float(err))
    else:
        for i, (r, y) in enumerate(zip(runs, vals)):
            x = pos_map[r] + 0.5
            g.SetPoint(i, float(x), float(y))
    g.SetMarkerStyle(marker_style)
    g.SetMarkerSize(marker_size)
    g.SetMarkerColor(marker_color)
    g.Draw("P SAME")

    if title_lines:
        t0 = ROOT.TLatex()
        t0.SetNDC(True)
        y0 = 0.85
        for j, line in enumerate(title_lines):
            t0.SetTextFont(62 if j == 0 else 42)
            t0.SetTextSize(0.045 if j == 0 else 0.035)
            t0.DrawLatex(0.12, y0 - 0.045*j, line)
        titles.append(t0)

    ROOT.gPad.RedrawAxis()
    c.Update()

    # 输出
    os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
    c.SaveAs(output_pdf)
    if png_also:
        c.SaveAs(output_pdf.replace(".pdf", ".png"))
    if root_also:
        rfdir = os.path.dirname(root_also)
        if rfdir: os.makedirs(rfdir, exist_ok=True)
        rf = ROOT.TFile.Open(root_also, "RECREATE")
        c.Write("canvas")
        rf.Close()

    # 返回对象，避免被 GC
    return {
        "canvas": c,
        "frame": frame,
        "graph": g,
        "boxes": boxes,
        "texts": texts,
        "titles": titles,
        "runs_seen": runs_seen
    }

# batch mode
ROOT.gROOT.SetBatch(True)

folder_301_waveform = "dump/300_RootConverterX/beamtests"
list_301_waveform_root_files = [f for f in os.listdir(folder_301_waveform) if f.endswith(".root") and not f.startswith("._")]
output_folder = "dump/QA_Results"
output_file_path = os.path.join(output_folder, "QA_301_Waveform_Results.root")

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
# file name: dump/301_Waveform/beamtests/Run0179.root
list_301_waveform_run_numbers = []
for filename in list_301_waveform_root_files:
    run_number_str = filename.replace("Run", "").replace(".root", "")
    try:
        run_number = int(run_number_str)
        list_301_waveform_run_numbers.append(run_number)
    except ValueError:
        continue
print("Found 301_Waveform root files for run numbers:", list_301_waveform_run_numbers)

# list_301_waveform_hamming_code_error_rates = []
# list_301_waveform_daqh_header_footer_error_rates = []
# list_301_waveform_toa_valid_fractions = []
list_301_waveform_kColors = []
list_301_data_lines_per_trigger_averages    = []
list_301_data_lines_per_trigger_errors      = []
list_301_illegal_40_timestamp_averages      = []
list_301_illegal_40_timestamp_errors        = []
list_301_trigger_160byte_parse_failures     = []
list_301_trigger_160byte_parse_failures_err = []
list_301_trigger_idle_lines_averages        = []
list_301_trigger_idle_lines_errors          = []
list_301_trigger_header_ts_errors_averages  = []
list_301_trigger_header_ts_errors_errors    = []
list_301_trigger_skipped_lines_averages     = []
list_301_trigger_skipped_lines_errors       = []
list_301_samples_total = []
list_301_samples_valid = []


# go through each run number and extract QA values
for run_index in range(len(list_301_waveform_run_numbers)):
    run_number = list_301_waveform_run_numbers[run_index]
    root_file_path = os.path.join(folder_301_waveform, f"Run{run_number:04d}.root")
    input_root = ROOT.TFile.Open(root_file_path, "READ")
    if not input_root or input_root.IsZombie():
        print(f"Failed to open root file: {root_file_path}")
        continue

    Avg_Data_Lines_Per_Trigger = None
    Err_Data_Lines_Per_Trigger = None
    Avg_Trigger_Illegal_40_Timestamps = None
    Err_Trigger_Illegal_40_Timestamps = None
    Avg_Trigger_160Byte_Parse_Failures = None
    Err_Trigger_160Byte_Parse_Failures = None
    Avg_Trigger_Idle_Lines = None
    Err_Trigger_Idle_Lines = None
    Avg_Trigger_Header_Timestamp_Errors = None
    Err_Trigger_Header_Timestamp_Errors = None
    Avg_Trigger_Skipped_Lines = None
    Err_Trigger_Skipped_Lines = None
    Total_Samples_In_Total = None
    Total_Valid_Samples = None


    param_avg_data_lines = input_root.Get("Avg_Data_Lines_Per_Trigger")
    if param_avg_data_lines:
        Avg_Data_Lines_Per_Trigger = param_avg_data_lines.GetVal()
    param_err_data_lines = input_root.Get("Err_Data_Lines_Per_Trigger")
    if param_err_data_lines:
        Err_Data_Lines_Per_Trigger = param_err_data_lines.GetVal()
    param_avg_illegal_40 = input_root.Get("Avg_Trigger_Illegal_40_Timestamps")
    if param_avg_illegal_40:
        Avg_Trigger_Illegal_40_Timestamps = param_avg_illegal_40.GetVal()
    param_err_illegal_40 = input_root.Get("Err_Trigger_Illegal_40_Timestamps")
    if param_err_illegal_40:
        Err_Trigger_Illegal_40_Timestamps = param_err_illegal_40.GetVal()
    param_avg_160byte = input_root.Get("Avg_Trigger_160Byte_Parse_Failures")
    if param_avg_160byte:
        Avg_Trigger_160Byte_Parse_Failures = param_avg_160byte.GetVal()
    param_err_160byte = input_root.Get("Err_Trigger_160Byte_Parse_Failures")
    if param_err_160byte:
        Err_Trigger_160Byte_Parse_Failures = param_err_160byte.GetVal()
    param_avg_idle_lines = input_root.Get("Avg_Trigger_Idle_Lines")
    if param_avg_idle_lines:
        Avg_Trigger_Idle_Lines = param_avg_idle_lines.GetVal()
    param_err_idle_lines = input_root.Get("Err_Trigger_Idle_Lines")
    if param_err_idle_lines:
        Err_Trigger_Idle_Lines = param_err_idle_lines.GetVal()
    param_avg_header_ts_errors = input_root.Get("Avg_Trigger_Header_Timestamp_Errors")
    if param_avg_header_ts_errors:
        Avg_Trigger_Header_Timestamp_Errors = param_avg_header_ts_errors.GetVal()
    param_err_header_ts_errors = input_root.Get("Err_Trigger_Header_Timestamp_Errors")
    if param_err_header_ts_errors:
        Err_Trigger_Header_Timestamp_Errors = param_err_header_ts_errors.GetVal()
    param_avg_skipped_lines = input_root.Get("Avg_Trigger_Skipped_Lines")
    if param_avg_skipped_lines:
        Avg_Trigger_Skipped_Lines = param_avg_skipped_lines.GetVal()
    param_err_skipped_lines = input_root.Get("Err_Trigger_Skipped_Lines")
    if param_err_skipped_lines:
        Err_Trigger_Skipped_Lines = param_err_skipped_lines.GetVal()

    Total_Samples_In_Total = input_root.Get("Total_Samples_In_Total").GetVal()
    Total_Valid_Samples = input_root.Get("Total_Valid_Samples").GetVal()

    list_301_data_lines_per_trigger_averages.append(Avg_Data_Lines_Per_Trigger)
    list_301_data_lines_per_trigger_errors.append(Err_Data_Lines_Per_Trigger)
    list_301_illegal_40_timestamp_averages.append(Avg_Trigger_Illegal_40_Timestamps)
    list_301_illegal_40_timestamp_errors.append(Err_Trigger_Illegal_40_Timestamps)
    list_301_trigger_160byte_parse_failures.append(Avg_Trigger_160Byte_Parse_Failures)
    list_301_trigger_160byte_parse_failures_err.append(Err_Trigger_160Byte_Parse_Failures)
    list_301_trigger_idle_lines_averages.append(Avg_Trigger_Idle_Lines)
    list_301_trigger_idle_lines_errors.append(Err_Trigger_Idle_Lines)
    list_301_trigger_header_ts_errors_averages.append(Avg_Trigger_Header_Timestamp_Errors)
    list_301_trigger_header_ts_errors_errors.append(Err_Trigger_Header_Timestamp_Errors)
    list_301_trigger_skipped_lines_averages.append(Avg_Trigger_Skipped_Lines)
    list_301_trigger_skipped_lines_errors.append(Err_Trigger_Skipped_Lines)

    list_301_samples_total.append(Total_Samples_In_Total)
    list_301_samples_valid.append(Total_Valid_Samples)

    list_301_waveform_kColors.append(channel_color_map.get(run_number, {}).get("color", ROOT.kGray+1))
    input_root.Close()

sorted_indices = sorted(range(len(list_301_waveform_run_numbers)), key=lambda i: list_301_waveform_run_numbers[i])
list_301_waveform_run_numbers = [list_301_waveform_run_numbers[i] for i in sorted_indices]
list_301_data_lines_per_trigger_averages    = [list_301_data_lines_per_trigger_averages[i] for i in sorted_indices]
list_301_data_lines_per_trigger_errors      = [list_301_data_lines_per_trigger_errors[i] for i in sorted_indices]
list_301_illegal_40_timestamp_averages      = [list_301_illegal_40_timestamp_averages[i] for i in sorted_indices]
list_301_illegal_40_timestamp_errors        = [list_301_illegal_40_timestamp_errors[i] for i in sorted_indices]
list_301_trigger_160byte_parse_failures     = [list_301_trigger_160byte_parse_failures[i] for i in sorted_indices]
list_301_trigger_160byte_parse_failures_err = [list_301_trigger_160byte_parse_failures_err[i] for i in sorted_indices]
list_301_trigger_idle_lines_averages        = [list_301_trigger_idle_lines_averages[i] for i in sorted_indices]
list_301_trigger_idle_lines_errors          = [list_301_trigger_idle_lines_errors[i] for i in sorted_indices]
list_301_trigger_header_ts_errors_averages  = [list_301_trigger_header_ts_errors_averages[i] for i in sorted_indices]
list_301_trigger_header_ts_errors_errors    = [list_301_trigger_header_ts_errors_errors[i] for i in sorted_indices]
list_301_trigger_skipped_lines_averages     = [list_301_trigger_skipped_lines_averages[i] for i in sorted_indices]
list_301_trigger_skipped_lines_errors       = [list_301_trigger_skipped_lines_errors[i] for i in sorted_indices]
list_301_samples_total = [list_301_samples_total[i] for i in sorted_indices]
list_301_samples_valid = [list_301_samples_valid[i] for i in sorted_indices]
list_301_waveform_kColors = [list_301_waveform_kColors[i] for i in sorted_indices]

channel_color_map = {}
for color_hex, run_info in color_json.items():
    run_list = run_info.get("run_numbers", [])
    r = int(color_hex[1:3], 16)
    g = int(color_hex[3:5], 16)
    b = int(color_hex[5:7], 16)
    color_code = ROOT.TColor.GetColor(r, g, b)
    description = run_info.get("description", "")
    for run_number in run_list:
        channel_color_map[run_number] = {
            "color": color_code,
            "description": description
        }

# 生成 run->color / run->text 两个 dict
run_to_color = {r: info["color"] for r, info in channel_color_map.items()}
run_to_text  = {r: info.get("description", "") for r, info in channel_color_map.items()}

# ========= 2) Draw the Data Lines per Trigger plot =========
title_lines = [
    "FoCal-H Prototype 3",
    "Beam Test 2025 Oct",
    "Data Lines per Trigger",
    time.strftime("Date: %Y-%m-%d")
]
out_pdf1 = os.path.join(output_folder, "QA_301_data_lines_per_trigger.pdf")

draw_run_categorical_bands_scatter(
    runs=list_301_waveform_run_numbers,
    y_values=list_301_data_lines_per_trigger_averages,
    y_errors=list_301_data_lines_per_trigger_errors,
    run_to_color=run_to_color,
    output_pdf=out_pdf1,
    y_min=0, y_max=max(list_301_data_lines_per_trigger_averages)*2.0, logy=False,
    band_shrink=0.00,
    y_label="Data Lines per Trigger",
    x_label="Run Number",
    title_lines=title_lines,
    run_to_text=run_to_text,
    png_also=True,
    root_also=os.path.join(output_folder, "QA_301_data_lines_per_trigger.root")
)

# ========= 3) Draw the Illegal 40 Timestamp per Trigger plot =========
title_lines = [
    "FoCal-H Prototype 3",
    "Beam Test 2025 Oct",
    "Illegal 40-byte Timestamp per Trigger",
    time.strftime("Date: %Y-%m-%d")
]
out_pdf2 = os.path.join(output_folder, "QA_301_illegal_40_timestamp_per_trigger.pdf")

draw_run_categorical_bands_scatter(
    runs=list_301_waveform_run_numbers,
    y_values=list_301_illegal_40_timestamp_averages,
    y_errors=list_301_illegal_40_timestamp_errors,
    run_to_color=run_to_color,
    output_pdf=out_pdf2,
    y_min=0, y_max=max(list_301_illegal_40_timestamp_averages)*2.0, logy=False,
    band_shrink=0.00,
    y_label="Illegal 40-byte Timestamp per Trigger",
    x_label="Run Number",
    title_lines=title_lines,
    run_to_text=run_to_text,
    png_also=True,
    root_also=os.path.join(output_folder, "QA_301_illegal_40_timestamp_per_trigger.root")
)

# ========= 4) Draw the 160-Byte Parse Failures per Trigger plot =========
title_lines = [
    "FoCal-H Prototype 3",
    "Beam Test 2025 Oct",
    "160-Byte Parse Failures per Trigger",
    time.strftime("Date: %Y-%m-%d")
]
out_pdf3 = os.path.join(output_folder, "QA_301_160byte_parse_failures_per_trigger.pdf")

draw_run_categorical_bands_scatter(
    runs=list_301_waveform_run_numbers,
    y_values=list_301_trigger_160byte_parse_failures,
    y_errors=list_301_trigger_160byte_parse_failures_err,
    run_to_color=run_to_color,
    output_pdf=out_pdf3,
    y_min=0, y_max=max(list_301_trigger_160byte_parse_failures)*2.0, logy=False,
    band_shrink=0.00,
    y_label="160-Byte Parse Failures per Trigger",
    x_label="Run Number",
    title_lines=title_lines,
    run_to_text=run_to_text,
    png_also=True,
    root_also=os.path.join(output_folder, "QA_301_160byte_parse_failures_per_trigger.root")
)   

# ========= 5) Draw the Idle Lines per Trigger plot =========
title_lines = [
    "FoCal-H Prototype 3",
    "Beam Test 2025 Oct",
    "Idle Lines per Trigger",
    time.strftime("Date: %Y-%m-%d")
]
out_pdf4 = os.path.join(output_folder, "QA_301_idle_lines_per_trigger.pdf")

y_max_idle = max([v for v in list_301_trigger_idle_lines_averages if v is not None] + [1.0]) * 2.0

draw_run_categorical_bands_scatter(
    runs=list_301_waveform_run_numbers,
    y_values=list_301_trigger_idle_lines_averages,
    y_errors=list_301_trigger_idle_lines_errors,
    run_to_color=run_to_color,
    output_pdf=out_pdf4,
    y_min=0, y_max=y_max_idle, logy=False,
    band_shrink=0.00,
    y_label="Idle Lines per Trigger",
    x_label="Run Number",
    title_lines=title_lines,
    run_to_text=run_to_text,
    png_also=True,
    root_also=os.path.join(output_folder, "QA_301_idle_lines_per_trigger.root")
)

# ========= 6) Draw the Header Timestamp Errors per Trigger plot =========
title_lines = [
    "FoCal-H Prototype 3",
    "Beam Test 2025 Oct",
    "Header Timestamp Errors per Trigger",
    time.strftime("Date: %Y-%m-%d")
]
out_pdf5 = os.path.join(output_folder, "QA_301_header_timestamp_errors_per_trigger.pdf")

draw_run_categorical_bands_scatter(
    runs=list_301_waveform_run_numbers,
    y_values=list_301_trigger_header_ts_errors_averages,
    y_errors=list_301_trigger_header_ts_errors_errors,
    run_to_color=run_to_color,
    output_pdf=out_pdf5,
    y_min=0, y_max=max([v for v in list_301_trigger_header_ts_errors_averages if v is not None] + [1.0]) * 2.0, logy=False,
    band_shrink=0.00,
    y_label="Header Timestamp Errors per Trigger",
    x_label="Run Number",
    title_lines=title_lines,
    run_to_text=run_to_text,
    png_also=True,
    root_also=os.path.join(output_folder, "QA_301_header_timestamp_errors_per_trigger.root")
)   

# ========= 7) Draw the Skipped Lines per Trigger plot =========
title_lines = [
    "FoCal-H Prototype 3",
    "Beam Test 2025 Oct",
    "Skipped Lines per Trigger",
    time.strftime("Date: %Y-%m-%d")
]
out_pdf6 = os.path.join(output_folder, "QA_301_skipped_lines_per_trigger.pdf")

draw_run_categorical_bands_scatter(
    runs=list_301_waveform_run_numbers,
    y_values=list_301_trigger_skipped_lines_averages,
    y_errors=list_301_trigger_skipped_lines_errors,
    run_to_color=run_to_color,
    output_pdf=out_pdf6,
    y_min=0, y_max=max([v for v in list_301_trigger_skipped_lines_averages if v is not None] + [1.0]) * 2.0, logy=False,
    band_shrink=0.00,
    y_label="Skipped Lines per Trigger",
    x_label="Run Number",
    title_lines=title_lines,
    run_to_text=run_to_text,
    png_also=True,
    root_also=os.path.join(output_folder, "QA_301_skipped_lines_per_trigger.root")
)

# ========= 8) Draw the Sample Numbers plot =========
title_lines = [
    "FoCal-H Prototype 3",
    "Beam Test 2025 Oct",
    "Total and Valid Sample Counts",
    time.strftime("Date: %Y-%m-%d")
]
out_pdf7 = os.path.join(output_folder, "QA_301_sample_counts.pdf")

y_lists = [list_301_samples_total, list_301_samples_valid]
y_error_lists = [[0]*len(list_301_samples_total), [0]*len(list_301_samples_valid)]

draw_run_categorical_bands_multi(
    runs=list_301_waveform_run_numbers,
    y_value_lists=y_lists,
    y_error_lists=y_error_lists,
    run_to_color=run_to_color,
    output_pdf=out_pdf7,
    y_min=0, y_max=max(list_301_samples_total)*1.2, logy=False,
    band_shrink=0.00,
    y_label="Sample Counts",
    x_label="Run Number",
    title_lines=title_lines,
    run_to_text=run_to_text,
    png_also=True,
    root_also=os.path.join(output_folder, "QA_301_sample_counts.root"),
    marker_colors=[ROOT.kBlue, ROOT.kGreen+2]
)