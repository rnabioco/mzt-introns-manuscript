wget https://raw.githubusercontent.com/kriemo/toolbucket/master/ena_downloader.py

python \
    ./ena_downloader.py \
    -s PRJNA342320 \
    -d ascp \
    -p ~/.aspera/connect/bin/ascp \
    -k ~/.aspera/connect/etc/asperaweb_id_dsa.openssh

