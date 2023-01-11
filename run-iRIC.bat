SET ROPTS=--no-save --no-environ --no-init-file --no-restore --no-Rconsole

start http://127.0.0.1:7777 & R-Portable\App\R-Portable\bin\x64\Rscript.exe iRIC.R 1> log_iRIC.log 2>&1