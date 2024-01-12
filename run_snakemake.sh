CLUSTER_CMD=("bsub -n {threads} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -e {cluster.error} -J {cluster.name} -W {cluster.time}")
config_yaml=metadata/configs.yaml
cluster_yaml=metadata/cluster.yaml

cmd="snakemake"
cmd="$cmd --configfile $config_yaml"
cmd="$cmd --jobs 1000"
cmd="$cmd --restart-times 0"
cmd="$cmd --rerun-incomplete"
cmd="$cmd --cluster-config $cluster_yaml"
cmd="$cmd --cluster \"${CLUSTER_CMD}\""
cmd="$cmd --cluster-cancel bkill"
cmd="$cmd --use-singularity"
cmd="$cmd -p"
cmd="$cmd --singularity-args \"--bind /juno --bind /rtsess01\""
cmd="$cmd --dry-run"

echo $cmd
eval $cmd
