<<<<<<< HEAD
=======
mkdir log
>>>>>>> 2bbac3071cc55c1d7853a1d1a92e810372a17376
echo "" > log/info.log
echo "" > log/err.log
timestamp=$(date +"%Y%m%d_%H%M%S")
nohup python src/main.py > log/info.log 2> log/err.log & 
pid=$!
echo $pid > log/pipeline_$timestamp.pid
echo "Pipeline started with PID: $pid at $timestamp"
