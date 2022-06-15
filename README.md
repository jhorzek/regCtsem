# regCtsem

Regularized Continuous Time Structural Equation Modeling. regCtsem builds on the ctsemOMX and ctsem package (Driver et al., 2017) and extends it with
different regularization techniques based on Jacobucci et al. (2016) and Huang et al. (2017).

# Installation

If you want to install regCtsem from GitHub, use the following commands in R:

    if(!require(devtools))install.packages("devtools")

    devtools::install_github("jhorzek/regCtsem")
    
In case you also want to install the vignettes, use:


    if(!require(devtools))install.packages("devtools")

    devtools::install_github("jhorzek/regCtsem", build_vignettes = TRUE)
    
    
# Getting Started

A good place to start is the help page of the main function:

    ?regCtsem::regCtsem

If you also installed the vignettes, use ``vignette("regCtsem", package = "regCtsem")'' to get some more examples and general guidelines.

# References

* Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modeling with R Package ctsem. Journal of Statistical Software, 77(5). https://doi.org/10.18637/jss.v077.i05
* Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
* Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. Structural Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793


# Important Notes

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 