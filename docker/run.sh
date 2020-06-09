docker run \
  -it \
  --rm \
  --name 16S_container \
  --mount type=bind,source=/Users/jordycoolen/PycharmProjects/16S,target=/16S \
  --mount type=bind,source=/Volumes/ELEMENTS/NimaGen/kraken2_db,target=/16S/kraken2_db/,readonly \
  --mount type=bind,source=/Volumes/ELEMENTS/centrifuge_db,target=/16S/centrifuge_db/,readonly \
  --mount type=bind,source=/Volumes/ELEMENTS/NimaGen/full_db/,target=/16S/full_db/,readonly \
  --mount type=bind,source=/Volumes/ELEMENTS/NimaGen/data,target=/seq,readonly \
  16s:latest /bin/bash