# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.


## Scientific context for the project

The paper is in https://arxiv.org/abs/2606.12597 . A (not committed, since useless) symlink to the latex file is provided too.


## R conventions for all my R work (this project or otherwise)

- Stick to standard R as much as possible and minimize use of the tidyverse. Specifically:

- Do not use dplyr. If data.frame is slow and we need an alternative, use data.table.

- Never use pipes (from tidyverse or the new ones in R).

- Do not use the stringr package. Use stringi.

- ggplot2 is OK, but respecting the above rules.

- Avoid (unless essential) any non-standard evaluation. We want clean, standard, non-surprises R code that will run 5 years from now.

- Do not use Roxygen2 or similar for documentation (you might see remnants are commented out). I write Rd files directly (not in this project) and use comments in the code.

- Never put 2 or more statements per line separated by a ";".

- Use spaces around "+", "-", "=", "*", etc.

- Except for one-line functions, always use "return" explicitly (i.e., do not just use the 'last thing evaluated before exiting is the return value').

- Never, ever, use global variables (except under extreme, very hard to handle otherwise, scenarios). In a file with functions called from other files, global variables are unacceptable — this covers both defining a global and depending on one: a function must get everything it needs from its arguments, never from a free variable expected to exist in the caller's environment. (In an interactive file, say to produce some figures, a global variable might be innocuous, but we don't have any of these here).

- Comments and line width: use the 80-column rule, but for comments especially, only up to column 74; so start a newline if you go beyond column 74. See my code for examples.

- In comments: do not use ALLCAPS for words. Do not write NEVER, but never (or Never), unless something can really kill you or similar.

- I use `checkUsageEnv(env = .GlobalEnv)` a lot, at the end of most files, to catch problems early.

### Comments
- Follow my standards: if possible (unless trivial) a function signature a la "How to design programs" and a minimal comment of what/why is done, right before the function. Details, inside the function. Rationale is: a user should get the details required to run the function, but the logic/tricks/shortcuts are something one reads only occasionally. There might be exceptions to this.
- Use clear, simple language in comments.

## Launching R

- When launching R so it reads a script I prefer this invocation:
  - `R --vanilla -f <input_file_name> &> <output_file_name>`.
  - That invocation: avoids all the `/bin/sh BATCH` layer, is clear that it sends both stdout and stderr to `output_file_name`, and overwrites the file (no appending).
- All these scripts should start with
```
date()
version
```
and finish with `date()`. If you use OncoSimulR or evamtools, also add, right after loading them `packageVersion("OncoSimulR")` or `packageVersion("evamtools")`.

  - For further redundancy, add `sessionInfo()`.

- Unless there is a strong reason to do otherwise, always run `input_file.R` from the directory where `input_file.R` lives.

- For intermediate tests you run (e.g., testing), you can use `Rscript` if you prefer. For real work that will be saved, use the above invocation. If you think the above invocation is a bad idea on a particular situation, discuss it with me.

- Write inline R code to a temp file instead of using Rscript -e '...' with multiline strings.

- If we need to use scripts that take arguments (processed with R's commandArgs) we will discuss how best to do it.


## Permissions

- Please, DO NOT ask for permission to read any file in this directory: You have my permission.

- Please, DO NOT ask for permission to execute R code that tests the files you create. You have my permission to run R and Rscript.

- Please, DO NOT ask for permission to run tests of R code. You have my permission.


## Multiple executions

- If we are running on the laptop, limit use of CPUs to 4 (`detectCores()` gives a bogus 14, but those aren't real floating-point decent CPUs, so things slow down); there is code here where number of cores is limited based on the hostname.

- On the laptop, never launch several concurrent executions, each of up to 4 cores.
