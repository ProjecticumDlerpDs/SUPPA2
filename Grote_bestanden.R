# Installeer Git LFS (als het nog niet geïnstalleerd is) via de terminal van macbook
brew install git-lfs
git lfs install

# Vertel Git welke bestanden met LFS moeten worden beheerd
git lfs track "*.csv.gz"
git lfs track "*.mtx.gz"

# Voeg de wijziging toe in .gitattributes
git add .gitattributes

# Commit en push
git commit -m "Add large count matrix via Git LFS"
git push origin main
