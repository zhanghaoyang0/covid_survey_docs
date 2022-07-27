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
git commit -m "20220722"
git remote -v 
# git remote rm origin # delelte origin, if have https link need password each time
# git remote add origin git@gitee.com:zhanghaoyang0/slurm_gwas_code.git 
git push -u origin "master"

git rev-list --objects --all | grep -E `git verify-pack -v .git/objects/pack/*.idx | sort -k 3 -n | tail -10 | awk '{print$1}' | sed ':a;N;$!ba;s/\n/|/g'`
git filter-branch --tree-filter 'rm -f /home/yanglab_data/user/zhanghy/project/slurm_gwas_codenafld_code/' --tag-name-filter cat -- --all
git push origin --tags --force
git push origin --all --force

git filter-branch --force --index-filter 'git rm --cached -r --ignore-unmatch nafld_code/E11_0.8_f.tsv' --prune-empty --tag-name-filter cat -- --all



# git pull origin "master"
### have repo
cd existing_git_repo
git remote add origin https://gitee.com/zhanghaoyang0/nafld.git
git push -u origin "master"
