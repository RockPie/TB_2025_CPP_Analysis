    # read_sheet.py
import os
from typing import Optional
import gspread
import pandas as pd
from google.oauth2.service_account import Credentials
import json
import paramiko
import stat
import posixpath

SHEET_TITLE = "Logbook_FoCalH_SPS_Oct2025"
WORKSHEET   = "Logbook"
SA_KEY      = "config/tb-2025-oct-automation-53c24706c7ea.json"
VERSION     = "0.1"

SCOPES = [
    "https://www.googleapis.com/auth/spreadsheets",   # <-- CHANGED (write)
    "https://www.googleapis.com/auth/drive.readonly",
]

def get_client(sa_key_path: Optional[str] = None) -> gspread.Client:
    """基于 Service Account 创建 gspread 客户端。"""
    if sa_key_path and os.path.exists(sa_key_path):
        creds = Credentials.from_service_account_file(sa_key_path, scopes=SCOPES)
    else:
        creds = Credentials.from_service_account_file(
            os.environ["GOOGLE_APPLICATION_CREDENTIALS"], scopes=SCOPES
        )
    return gspread.authorize(creds)

def open_ws(client: gspread.Client, sheet_title: str, worksheet_name: str):
    sh = client.open(sheet_title)
    return sh.worksheet(worksheet_name)

def get_header_map(ws) -> dict:
    header = ws.row_values(1)
    return {name.strip(): idx+1 for idx, name in enumerate(header)}

def find_row_by_run(ws, run_number: int, run_col_name="Run Number") -> int:
    hdr = get_header_map(ws)
    if run_col_name not in hdr:
        return -1
    col_idx = hdr[run_col_name]
    col_vals = ws.col_values(col_idx)[1:]
    for i, v in enumerate(col_vals, start=2):
        try:
            if int(str(v).strip()) == int(run_number):
                return i
        except ValueError:
            continue
    return -1

def update_sheet_after_upload(ws, run_number: int, remote_root_path: str, version: str):
    """把 Root File Path / Root File Version 写回到对应行。"""
    hdr = get_header_map(ws)
    row = find_row_by_run(ws, run_number)
    if row == -1:
        print(f"[WARN] Run Number={run_number} not found in sheet; skip writing back.")
        return

    path_col_candidates = ["Root File Path", "RootFilePath"]
    ver_col_candidates  = ["Root File Version", "RootFileVersion"]

    path_col = next((hdr[c] for c in path_col_candidates if c in hdr), None)
    ver_col  = next((hdr[c] for c in ver_col_candidates  if c in hdr), None)

    updates = []
    if path_col:
        updates.append({"range": gspread.utils.rowcol_to_a1(row, path_col), "values": [[remote_root_path]]})
    if ver_col and version is not None:
        updates.append({"range": gspread.utils.rowcol_to_a1(row, ver_col), "values": [[str(version)]]})

    if updates:
        ws.batch_update([{"range": u["range"], "values": u["values"]} for u in updates])
        print(f"[OK] Sheet updated for Run {run_number}: path+version")
    else:
        print(f"[WARN] No 'Root File Path'/'Root File Version' columns found; nothing written.")

def load_cred(path):
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    for k in ("user", "host", "password"):
        if k not in d or not d[k]:
            raise ValueError(f"Missing '{k}' in {path}")
    return d["user"], d["host"], d["password"]

def sftp_connect(user, host, password, port=22):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(hostname=host, port=port, username=user, password=password, timeout=10)
    return client, client.open_sftp()

def sftp_exists(sftp, path):
    try:
        sftp.stat(path)
        return True
    except FileNotFoundError:
        return False

def sftp_isdir(sftp, path):
    try:
        return stat.S_ISDIR(sftp.stat(path).st_mode)
    except FileNotFoundError:
        return False

def sftp_makedirs(sftp, remote_dir):
    parts = []
    head = remote_dir
    while True:
        head, tail = posixpath.split(head)
        if tail:
            parts.append(tail)
        else:
            if head:
                parts.append(head)
            break
    parts = list(reversed(parts))
    path = ""
    for p in parts:
        if p == "":
            continue
        if path == "" and remote_dir.startswith("/"):
            path = "/" + p if p != "/" else "/"
        else:
            path = posixpath.join(path, p)
        try:
            if not sftp_isdir(sftp, path):
                sftp.mkdir(path)
        except IOError as e:
            if not sftp_isdir(sftp, path):
                raise e

def sftp_remove_file_if_exists(sftp, remote_path):
    try:
        st = sftp.stat(remote_path)
        if stat.S_ISDIR(st.st_mode):
            raise IsADirectoryError(f"{remote_path} is a directory, not removing.")
        sftp.remove(remote_path)
        print(f"[INFO] Removed existing file: {remote_path}")
    except FileNotFoundError:
        pass

def human(n):
    for unit in ["B","KB","MB","GB","TB"]:
        if n < 1024 or unit == "TB":
            return f"{n:.2f} {unit}"
        n /= 1024

def progress_callback(filename, bytes_transferred, bytes_total):
    pct = (bytes_transferred / bytes_total * 100) if bytes_total else 0
    end = "\n" if bytes_transferred == bytes_total else "\r"
    print(f"[UPLOAD] {filename}: {human(bytes_transferred)}/{human(bytes_total)} ({pct:5.1f}%)", end=end, flush=True)

def read_worksheet_as_dataframe(client: gspread.Client, sheet_title: str, worksheet_name: str) -> pd.DataFrame:
    sh = client.open(sheet_title)
    ws = sh.worksheet(worksheet_name)
    rows = ws.get_all_records(numericise_ignore=['all'])
    df = pd.DataFrame(rows)

    for col in ["Run", "Run Number", "RunNumber", "run", "run_number"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    for col in ["Good", "good", "IsGood"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.upper().map({'GOOD': 'GOOD', 'BAD': 'BAD'}).fillna('UNKNOWN')

    # for remote path
    for col in ["Remote Path", "remote_path", "RemotePath"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()

    # for 'Root File Path', 'Root File Version'
    for col in ["Root File Path", "RootFilePath", "Root File Version", "RootFileVersion"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()

    return df

def do_root_conversion(run_number: int, remote_path: str):
    import subprocess

    target_path = f"./dump/102_EventMatch/beamtests/Run{run_number:04d}.root"
    cmd = ["snakemake", "--cores", "4", "--rerun-incomplete", target_path]         # <-- CHANGED
    print("Executing:", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Conversion failed for Run {run_number}")
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        return None

    print(f"Conversion successful for Run {run_number}")

    if not os.path.exists(target_path):
        print(f"[ERROR] Expected output does not exist: {target_path}")
        return None


    remote_path_beam_dir = posixpath.dirname(remote_path)
    remote_path_data_dir = posixpath.dirname(remote_path_beam_dir)
    remote_path_focalh_dir = posixpath.dirname(remote_path_data_dir)
    remote_root_dir = posixpath.join(remote_path_focalh_dir, "root", "beam")  # <-- CHANGED to posixpath
    remote_file_name = f"Run{run_number:04d}.root"
    remote_file_path = posixpath.join(remote_root_dir, remote_file_name)       # <-- CHANGED to posixpath

    print(f"Uploading to {remote_file_path} ...")
    lxplus_cred = "config/lxplus.json"
    user, host, password = load_cred(lxplus_cred)
    ssh_client, sftp = sftp_connect(user, host, password)
    try:
        sftp_makedirs(sftp, remote_root_dir)
        sftp_remove_file_if_exists(sftp, remote_file_path)
        local_file_size = os.path.getsize(target_path)
        sftp.put(target_path, remote_file_path, callback=lambda x,y: progress_callback(remote_file_name, x, y))
        print(f"Upload completed: {remote_file_path} ({human(local_file_size)})")
    finally:
        sftp.close()
        ssh_client.close()
        # check if there is downloaded raw data in data/beamtests, if yes, remove it to save space
        # raw_data_path = f"./data/beamtests/Run{run_number:04d}.ch2g"
        # if os.path.exists(raw_data_path):
        #     os.remove(raw_data_path)
        #     print(f"Removed local raw data file: {raw_data_path}")

    return remote_file_path  # <-- ensure we return it on success
    # result = subprocess.run(cmd, capture_output=True, text=True)
    # if result.returncode == 0:
    #     print(f"Conversion successful for Run {run_number}")
    # else:
    #     print(f"Conversion failed for Run {run_number}")
    #     print("STDOUT:", result.stdout)
    #     print("STDERR:", result.stderr)

def main():
    client = get_client(SA_KEY if SA_KEY else None)
    ws = open_ws(client, SHEET_TITLE, WORKSHEET)               # <-- NEW
    df = read_worksheet_as_dataframe(client, SHEET_TITLE, WORKSHEET)

    print("Columns:", list(df.columns))
    print("Total rows:", len(df))

    if "Run Number" in df.columns and "Good" in df.columns:
        good_runs = df[df["Good"] == 'GOOD']
        for _, row in good_runs.iterrows():
            enable_conversion = False
            this_script_version = float(VERSION) if VERSION else 0.0
            row_root_path = str(row.get("Root File Path", "") or "").strip()
            row_root_version = str(row.get("Root File Version", "") or "").strip()
            print(f"Run {row['Run Number']}: Root Path='{row_root_path}', Version='{row_root_version}'")

            if (row_root_path == "") or (row_root_version and float(row_root_version) < this_script_version):
                enable_conversion = True

            if not enable_conversion:
                print(f"Skipping Run {row['Run Number']}: Already up-to-date.")
                continue

            run_number = int(row["Run Number"])
            remote_path = str(row.get("Remote Path", "") or "").strip()
            if not remote_path:
                print(f"Skipping Run {run_number}: No remote path provided.")
                continue

            remote_root_path = do_root_conversion(run_number, remote_path)
            if remote_root_path:
                # === 写回 Sheet：Root File Path / Root File Version ===
                update_sheet_after_upload(ws, run_number, remote_root_path, VERSION)   # <-- NEW

if __name__ == "__main__":
    main()