## global setting
git config --global user.name "zhanghy"
git config --global user.email "10880429+zhanghaoyang0@user.noreply.gitee.com"

### creat repo
mkdir slurm_gwas_code
cd slurm_gwas_code
git init 
touch README.md
git add README.md
git add .
git commit -m "20220525"
# git remote add origin git@gitee.com:zhanghaoyang0/slurm_gwas_code.git # if have orgin, use git remote rm origin to delelte. | if have https, need password each time
# git remote -v 
git push -u origin "master"

### have repo
cd existing_git_repo
git remote add origin https://gitee.com/zhanghaoyang0/nafld.git
git push -u origin "master"
