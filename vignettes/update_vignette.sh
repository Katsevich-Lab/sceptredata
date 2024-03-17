cp ~/research_code/sceptre-book/sceptredata.qmd ~/research_code/sceptredata/vignettes/sceptredata.Rmd

sed -i -e '/::: callout-note/,/:::/d' sceptredata.Rmd
sed -i -e 's/# Example data in the `sceptredata` package {#sec-sceptre_data}//g' sceptredata.Rmd
sed -i -e 's|\[@gasperini2019\]|\[(Gasperini, 2019)\](https://pubmed.ncbi.nlm.nih.gov/30612741/)|g' sceptredata.Rmd
sed -i -e 's|\[@papalexi2021\]|\[(Papalexi, 2021)\](https://pubmed.ncbi.nlm.nih.gov/33649593/)|g' sceptredata.Rmd
sed -i -e 's/{#[^}]*}//g' sceptredata.Rmd
sed -i -e 's/[^.]*@sec-[^.]*\.//g' sceptredata.Rmd
echo "
\`\`\`{r}
library(sessioninfo); session_info()
\`\`\`" >> sceptredata.Rmd
sed -i '' '1i\
--- \
title: "Example data in the `sceptredata` package" \
output: rmarkdown::html_vignette \
vignette: > \
  %\\VignetteIndexEntry{Example data in the sceptredata package} \
  %\\VignetteEngine{knitr::rmarkdown} \
  %\\VignetteEncoding{UTF-8} \
--- \
' sceptredata.Rmd

rm sceptredata.Rmd-e
