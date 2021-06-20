# pirouette_example_29

Branch   |[![Travis CI logo](pics/TravisCI.png)](https://travis-ci.com)                                                                                                 |[![AppVeyor logo](pics/AppVeyor.png)](https://appveyor.com)                                                                                               
---------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------
`master` |[![Build Status](https://travis-ci.com/richelbilderbeek/pirouette_example_29.svg?branch=master)](https://travis-ci.com/richelbilderbeek/pirouette_example_29) |[![Build status](https://ci.appveyor.com/api/projects/status/uh6ek769fcap0ydr/branch/master?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/pirouette-example-29/branch/master)
`develop`|[![Build Status](https://travis-ci.com/richelbilderbeek/pirouette_example_29.svg?branch=develop)](https://travis-ci.com/richelbilderbeek/pirouette_example_29)|[![Build status](https://ci.appveyor.com/api/projects/status/uh6ek769fcap0ydr/branch/develop?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/pirouette-example-29/branch/develop)

A [pirouette example](https://github.com/richelbilderbeek/pirouette_examples)
that shows the effect of MCMC chain length on ESS.

## Running on Peregrine

Install `pirouette` using the [peregrine](https://github.com/richelbilderbeek/peregrine)
bash and R scripts.

Then, in the main folder of this repo, type:

```
sbatch scripts/rerun.sh
```

## Results

 * Download the intermediate data at 
   [https://www.richelbilderbeek.nl/pirouette_example_29.zip](https://www.richelbilderbeek.nl/pirouette_example_29.zip)

![](errors_1.png)

![](errors_2.png)

![](errors_3.png)

![](errors_4.png)

