with open("/uufs/chpc.utah.edu/common/home/saarman-group1/uphlfiles/SPAdes_scripts/blast_result/one-hit-test-bold-chordata.tsv", "r+", encoding="utf-8", errors="ignore") as f:
    data = f.read().replace("|", "\t")
    f.seek(0)
    f.write(data)
    f.truncate()
