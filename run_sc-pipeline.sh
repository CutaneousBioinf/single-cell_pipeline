timestamp=$(date +"%Y%m%d_%H%M%S")
nohup python src/main.py > log/info.log 2> log/err.log & 
pid=$!
echo $pid > log/pipeline_$timestamp.pid
echo "Pipeline started with PID: $pid at $timestamp"