import subprocess
import google.datalab.storage as storage # pip install datalab
import os

uganda2k = [o.key for o in storage.Bucket('african-seq-data').objects() if o.key.startswith('uganda2k')]

with open('inputSamplesFile.txt', 'w') as f:
    for x in uganda2k:
        if x.endswith(".bam") or x.endswith(".cram"):
            file_path = "gs://african-seq-data/{}".format(x)
            base = os.path.basename(file_path) # get file name without path
            # get file name without .{bam,cram} extension: for ubam file
            file_no_ext = os.path.splitext(base)[0]
            file_no_ext = file_no_ext.replace("#", "_") # replace the # with _ because MergeVCFs in picard has a trouble dealing with special characters like #
            # e.g before this change, one of the file names was EGAR00001140731_10256_1#73 and will be changed to EGAR00001140731_10256_1_73
            sample = subprocess.run(("gsutil cat {} | samtools view -H | grep ^@RG | tr '\t' '\n' | grep -m1 '^SM:' | cut -d ':' -f 2").format(file_path),
                                    shell=True, capture_output=True, encoding="utf-8").stdout
            f.write(file_path + "\t" + file_no_ext + "\t" + sample + "\n")


