
Version 1.1.6 ----------------------------------------------------------
Dear reviewer,

This version removes `superheat` as a dependency.

Best,
Salvador

Version 1.1.5 ----------------------------------------------------------

Dear reviewer,

This version addresses the issues that might occur due to changes in ggplot2.

Version 1.1.4 ----------------------------------------------------------

## First resubmission.

R: Dear reviewer,

The raised issue has been addressed.

Best,
Salvador.

------------------------------------------------------------------------

Q: Thanks, we see:

   The Title field starts with the package name.

and "A Package" is redundant, so please write this as

Title: Complementary Genomic Functions

Please fix and resubmit.

Best,
Uwe Ligges

-----------------------------------------------------------------------

## Original submission

Dear reviewer,

We appreciate your time to verify our library.
It has been reviewed and re-tested. The check results are shown below.
Also, please find bellow the answers to your concerns and suggestions.

Best,

Salvador

-----------------------------------------------------------------------

Q: Dear maintainer,

Please see the problems shown on
<https://cran.r-project.org/web/checks/check_results_ASRgenomics.html>.

Specifically, pls see the Rd \usage sections check NOTEs.

Please correct before 2024-02-19 to safely retain your package on CRAN.

Best,
-k

-----------------------------------------------------------------------

## R CMD check results

## Test environments 

# - Windows Server 2022, R-devel, 64 bit

❯ checking CRAN incoming feasibility ... [14s] NOTE
❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''
❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'
    
0 errors ✔ | 0 warnings ✔ | 3 notes ✖

# - R-hub Ubuntu Linux 20.04.1 LTS, R-release, GCC

Build ID:   ASRgenomics_1.1.4.tar.gz-96b148dd6c584b05a9ef102b6692be06
Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC
Submitted:  12h 22m 59.7s ago
Build time: 3h 6m 39.4s

Build failed during preparation or aborted

# - R-hub Fedora Linux, R-devel, clang, gfortran

Build ID:   ASRgenomics_1.1.4.tar.gz-4d9c3cf847e944349cbf059c7280cc6d
Platform:   Fedora Linux, R-devel, clang, gfortran
Submitted:  12h 22m 59.7s ago
Build time: 2h 52m 13.3s

Build failed during preparation or aborted

# - winbuilder (devel)

* checking CRAN incoming feasibility ... [12s] NOTE

Version 1.1.3 ----------------------------------------------------------

Dear reviewer,

We appreciate your time to verify our library.
It has been reviewed and re-tested. The check results are shown below.
Also, please find bellow the answers to your concerns and suggestions.

Best,

Salvador

-----------------------------------------------------------------------

Q: Dear maintainer,

Please see the problems shown on
<https://cran.r-project.org/web/checks/check_results_ASRgenomics.html>.

Please correct before 2022-12-08 to safely retain your package on CRAN.

Do remember to look at the 'Additional issues'.

The CRAN Team

R: The failing tests were improved to avoid this issue. Also, one bug fix
has been done.

-----------------------------------------------------------------------

## R CMD check results

## Test environments 

# - R-hub windows-x86_64-devel (r-devel)

R-hub is failing with the error:

> * checking CRAN incoming feasibility ...
> Error in aspell(files, filter = list("dcf", ignore = ignore), control = control, :
> No suitable spell-checker program found

Hence, the library was locally tested under the following platform:

R version 4.2.2 (2022-11-25 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

And generates the following results:

0 errors | 0 warnings | 0 note

# - R-hub Ubuntu Linux 20.04.1 LTS, R-release, GCC

0 errors | 0 warnings | 1 note

* Days since last update: 3

# - R-hub Fedora Linux, R-devel, clang, gfortran

0 errors | 0 warnings | 1 note

* Days since last update: 3

# - winbuilder (devel)

0 errors | 0 warnings | 1 note

* Days since last update: 3

Version 1.1.2 ----------------------------------------------------------

## Second resubmission

Dear Victoria,

Thank you once again for your time and assistance.
The library has been reviewed according to the
suggested modifications.

Best,

Salvador

Version: 1.1.2
Date: 2022-11-17 18:21:22 UTC
SHA: 9765c5c391042f287b6a85b0bdeb82cac36a1d93

-----------------------------------------------------------------------

Q: I think a "the" is missing between "after" and "use" in "... before and 
after use of 'asreml' or another library to..."
and a "are" between "examples" and "decribed" in "... Association 
Analyses. Methods and examples described in Gezan, Oliveira, ..."

R: The suggested modifications have been applied.

-----------------------------------------------------------------------

## Resubmission

Dear Victoria,

We appreciate your time and effort to help us improve this library.
It has been reviewed and re-tested.
The check results are still the same from the initial submission.
Please find bellow the answers to your concerns and suggestions.

Best,

Salvador

Version: 1.1.2
Date: 2022-11-15 19:19:16 UTC
SHA: dc7384f8cae0ae68a8e0d5869383e05813c8492a

-----------------------------------------------------------------------

Q: If there are references describing the methods in your package, please 
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for 
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

R: This has been addressed by adding the manual to the referred field.

-----------------------------------------------------------------------

Q: "Using foo:::f instead of foo::f allows access to unexported objects. 
This is generally not recommended, as the semantics of unexported 
objects may be changed by the package author in routine maintenance."
Please omit one colon.
Used ::: in documentation:
      man/snp.pruning.Rd:
         map <- ASRgenomics:::dummy.map_(colnames(M.clean))

You have examples for unexported functions. Please either omit these 
examples or export these functions.
Examples for unexported function
   dummy.map_() in:
      snp.pruning.Rd

R: The section of the examples using internal functions have been removed.

-----------------------------------------------------------------------

Q: \dontrun{} should only be used if the example really cannot be executed 
(e.g. because of missing additional software, missing API keys, ...) by 
the user. That's why wrapping examples in \dontrun{} adds the comment 
("# Not run:") as a warning for the user.
Does not seem necessary.
Please unwrap the examples if they are executable in < 5 sec, or create 
additionally small toy examples to allow automatic testing.
(You could also replace \dontrun{} with \donttest, if it takes longer 
than 5 sec to be executed, but it would be preferable to have automatic 
checks for functions. Otherwise, you can also write some tests.)

R: The examples depending on external software have been maintained with 
\dontrun{}; the examples that run in less than five seconds have had
\dontrun{} removed; and the examples that take more than five seconds
have have \dontrun{} replaced with \donttest{}. We would like to keep
the examples as they are (not toy examples), because they tackle specific
'issues' found in the datasets and solve them.

-----------------------------------------------------------------------

Q: Please do not modify the global environment (e.g. by using <<-) in your 
functions. This is not allowed by the CRAN policies. e.g.: man/G.predict.Rd

R: The example has been modified.

-----------------------------------------------------------------------

## First submission

## Preamble

ASRgenomics has been available since 2020 via VSN International's website.
It is currently on version 1.1.2 which we wish to submit for consideration.

## R CMD check results

## Test environments 

# - R-hub windows-x86_64-devel (r-devel)

R-hub is failing with the error:

> * checking CRAN incoming feasibility ...
> Error in aspell(files, filter = list("dcf", ignore = ignore), control = control, :
> No suitable spell-checker program found

Hence, the library was locally tested under the following platform:

R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044

And generates the following results:

0 errors | 0 warnings | 1 note

* This is a new release.

# - R-hub Ubuntu Linux 20.04.1 LTS, R-release, GCC

0 errors | 0 warnings | 1 note

* This is a new release.

# - R-hub Fedora Linux, R-devel, clang, gfortran

0 errors | 0 warnings | 1 note

* This is a new release.
