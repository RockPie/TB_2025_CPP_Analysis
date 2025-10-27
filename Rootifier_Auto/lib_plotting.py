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
