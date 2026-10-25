```bash
cd /Users/golden/kongdd/R
Rscript -e 'install.packages("dplyr", repos = "https://cloud.r-project.org")'
Rscript final_sricpt/final_script.R
git status
git diff
git log --oneline -5
git add .
git commit -m "Update analysis"
git push
```
