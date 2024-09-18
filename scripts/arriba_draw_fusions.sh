set -veo pipefail
apptainer exec $@ docker://uhrigs/arriba:2.4.0 draw_fusions.sh