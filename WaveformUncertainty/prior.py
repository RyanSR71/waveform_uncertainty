import numpy as np
import bilby
import random
import time as tm
import sys
import scipy
import lal
from bilby.core import utils
from bilby.core.series import CoupledTimeAndFrequencySeries
from bilby.core.utils import PropertyAccessor
from bilby.gw.conversion import convert_to_lal_binary_neutron_star_parameters



def dphi_prior(phase_uncertainty,k, **kwargs):
    '''
    Generates a Gaussian prior for the phase correction parameters
    
    Parameters
    ===================
    phase_uncertainty: numpy.ndarray
        array of standard deviation of a set of phase differences; by default, this should be as a function of dimensionless frequency, xi
    k: int
        number of phase correction parameters desired
    mean_phase_difference: numpy.ndarray, optional
        array of the means of a set of phase differences, by default, this should be as a function of dimensionless frequency, xi
        default: None
    prior: bilby.core.prior.PriorDict, optional
        bilby prior object; if given, dphi priors will be added to this dictionary
        default: None
    geometrized: bool, optional
        if True, will return geometrized frequency nodes; if False, normal frequency nodes (Hz)
        default: True
    xi_low: float, optional
        if geometrized is True; lower bound on the geometrized frequency band
        default: 0.018
    xi_high: float, optional
        if geometrized is True; upper bound on the geometrized frequency band
        default: 1/pi (0.318...)
    f_low: float, optional
        if geometrized is False; lower bound on the normal frequency band
        default: 20.0 Hz
    f_high: float, optional
        if geometrized is False; upper bound on the normal frequency band
        default: 1024.0 Hz
        
    Returns
    ==================
    frequency_nodes: numpy.ndarray
        array of frequency nodes
    prior: bilby.core.prior.PriorDict
        bilby prior object containing the phase correction priors
    '''
    f_low = kwargs.get('f_low',20)
    f_high = kwargs.get('f_high',1024)
    xi_low = kwargs.get('xi_low',0.018)
    xi_high = kwargs.get('xi_high',1/np.pi)
    prior = kwargs.get('prior',None)
    geometrized = kwargs.get('geometrized',True)
    mean_phase_difference = kwargs.get('mean_phase_difference',None)
    
    if prior is None:
        prior = bilby.core.prior.PriorDict()
    
    if mean_phase_difference is None:
        mean_phase_difference = np.array([0]*len(phase_uncertainty))
    
    if geometrized is True:
        frequency_grid = np.linspace(0.001,1,len(phase_uncertainty))
        desired_frequency_nodes = np.geomspace(xi_low,xi_high,k+1)
    else:
        frequency_grid = np.linspace(f_low,f_high,len(phase_uncertainty))
        desired_frequency_nodes = np.geomspace(f_low,f_high,k+1)
        
    indexes = [list(frequency_grid).index(min(frequency_grid, key=lambda x:np.abs(x-node))) for node in desired_frequency_nodes]
    frequency_nodes = np.array(frequency_grid[indexes])

    prior['dphi_0'] = bilby.core.prior.DeltaFunction(name='dphi_0',latex_label=r'$\varphi_0$',peak=0)
    for i in list(range(len(frequency_nodes)))[1:]:
        prior[f'dphi_{i}'] = bilby.core.prior.Gaussian(name=f'dphi_{i}',latex_label=r'$\varphi_num$'.replace('num',str(i)),
                                                     mu=mean_phase_difference[indexes[i]],sigma=phase_uncertainty[indexes[i]])
    
    return frequency_nodes, prior



def dA_prior(amplitude_uncertainty,k, **kwargs):
    '''
    Generates a Gaussian prior for the amplitude correction parameters
    
    Parameters
    ===================
    amplitude_uncertainty: numpy.ndarray
        array of standard deviation of a set of amplitude differences; by default, this should be as a function of dimensionless frequency, xi
    k: int
        number of amplitude correction parameters desired
    mean_amplitude_difference: numpy.ndarray, optional
        array of the means of a set of amplitude differences, by default, this should be as a function of dimensionless frequency, xi
        default: None
    prior: bilby.core.prior.PriorDict, optional
        bilby prior object; if given, dA priors will be added to this dictionary
        default: None
    geometrized: bool, optional
        if True, will return geometrized frequency nodes; if False, normal frequency nodes (Hz)
        default: True
    xi_low: float, optional
        if geometrized is True; lower bound on the geometrized frequency band
        default: 0.018
    xi_high: float, optional
        if geometrized is True; upper bound on the geometrized frequency band
        default: 1/pi (0.318...)
    f_low: float, optional
        if geometrized is False; lower bound on the normal frequency band
        default: 20.0 Hz
    f_high: float, optional
        if geometrized is False; upper bound on the normal frequency band
        default: 1024.0 Hz
        
    Returns
    ==================
    frequency_nodes: numpy.ndarray
        array of frequency nodes
    prior: bilby.core.prior.PriorDict
        bilby prior object containing the amplitude correction priors
    '''
    f_low = kwargs.get('f_low',20)
    f_high = kwargs.get('f_high',1024)
    xi_low = kwargs.get('xi_low',0.018)
    xi_high = kwargs.get('xi_high',1/np.pi)
    prior = kwargs.get('prior',None)
    geometrized = kwargs.get('geometrized',True)
    mean_amplitude_difference = kwargs.get('mean_amplitude_difference',None)
    
    if prior is None:
        prior = bilby.core.prior.PriorDict()
    
    if mean_amplitude_difference is None:
        mean_amplitude_difference = np.array([0]*len(amplitude_uncertainty))
    
    if geometrized is True:
        frequency_grid = np.linspace(0.001,1,len(amplitude_uncertainty))
        desired_frequency_nodes = np.geomspace(xi_low,xi_high,k+1)
    else:
        frequency_grid = np.linspace(f_low,f_high,len(amplitude_uncertainty))
        desired_frequency_nodes = np.geomspace(f_low,f_high,k+1)
        
    indexes = [list(frequency_grid).index(min(frequency_grid, key=lambda x:np.abs(x-node))) for node in desired_frequency_nodes]
    frequency_nodes = np.array(frequency_grid[indexes])

    prior['dA_0'] = bilby.core.prior.DeltaFunction(name='dA_0',latex_label=r'$\alpha_0$',peak=0)
    for i in list(range(len(frequency_nodes)))[1:]:
        prior[f'dA_{i}'] = bilby.core.prior.Gaussian(name=f'dA_{i}',latex_label=r'$\alpha_num$'.replace('num',str(i)),
                                                     mu=mean_amplitude_difference[indexes[i]],sigma=amplitude_uncertainty[indexes[i]])
    
    return frequency_nodes, prior


def TotalMassConstraint(*,name,f_low,f_high,**kwargs):
    '''
    Generates a bilby prior to constrain the total mass
    
    Parameters
    ===================
    name: string
        name of the prior
    f_low: float
        lower bound on the frequency band
    f_high: float
        upper bound of the frequency band
    unit: string, optional
        unit of the parameter
        default: r'$\mathrm{M}_\odot$' (solar mass)
    latex_label: string, optional
        label for the parameter in LaTeX
        default: r'$M$'
    boundary: string, optional
        boundary condition type for the prior
        default: None
    xi_low: float, optional
        lower bound on the waveform uncertainty correction in dimensionless frequency
        default: 0.018
    xi_high: float, optional
        upper bound on the waveform uncertainty correction in dimensionless frequency
        default: 1/pi, 0.318...

    Returns
    ==================
    total_mass_prior: bilby.core.prior.base.Constraint
        bilby constraint prior object for the total mass
    '''
    unit = kwargs.get('unit',r'$\mathrm{M}_\odot$')
    latex_label = kwargs.get('latex_label',r'$M$')
    boundary = kwargs.get('boundary',None)
    xi_low = kwargs.get('xi_low',0.018)
    xi_high = kwargs.get('xi_high',1/np.pi)
    
    total_mass_prior = bilby.core.prior.Constraint(name=name,latex_label=latex_label,minimum=xi_high*203025.4467280836/f_high, maximum=xi_low*203025.4467280836/f_low, unit=unit)
    
    return total_mass_prior


def DeltaFConstraint(*,name,duration,f_low,f_high,**kwargs):
    '''
    Generates a bilby prior to constrain the total mass
    
    Parameters
    ===================
    name: string
        name of the prior
    duration: float
        duration of the signal
    f_low: float
        lower bound on the frequency band
    f_high: float
        upper bound of the frequency band
    n: int
        number of frequency nodes
    unit: string, optional
        unit of the parameter
        default: None
    latex_label: string, optional
        label for the parameter in LaTeX
        default: r'$\Delta f$'
    boundary: string, optional
        boundary condition type for the prior
        default: None

    Returns
    ==================
    total_mass_prior: bilby.core.prior.base.Constraint
        bilby constraint prior object for delta_f
    '''
    unit = kwargs.get('unit',None)
    latex_label = kwargs.get('latex_label',r'$\Delta f$')
    boundary = kwargs.get('boundary',None)
    
    delta_f_prior = bilby.core.prior.Constraint(name=name,latex_label=latex_label,minimum=4/duration, maximum=f_high-f_low, unit=unit)
    
    return delta_f_prior



def total_mass_conversion(parameters):
    '''
    Conversion function to generate the total mass from a set of parameters; to be used alongside the total mass prior
    
    Parameters
    ==================
    parameters: dict
        dictionary of binary black hole parameters
    
    Returns
    ==================
    parameters: dict
        input parameters, but with the total mass added
    '''
    parameters['total_mass'] = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
    
    return parameters


def delta_f_conversion(parameters,conversion_arguments):
    '''
    Conversion function to generate the total mass from a set of parameters; to be used alongside the total mass prior
    
    Parameters
    ==================
    parameters: dict
        dictionary of binary black hole parameters
    conversion_arguments: dict
        dictionary of additional arguments for the conversion function
    
    Returns
    ==================
    parameters: dict
        input parameters, but with delta_f added
    '''
        
    n = conversion_arguments['n']
    total_mass = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
    parameters['delta_f'] = (203025.4467280836/total_mass)*(parameters['xi_low']**(1-1/n)*parameters['xi_high']**(1/n)-parameters['xi_low'])
    return parameters

def bilby_delta_f_conversion(parameters):
    '''
    Conversion function to generate the total mass from a set of parameters; to be used alongside the total mass prior
    
    Parameters
    ==================
    parameters: dict
        dictionary of binary black hole parameters
    conversion_arguments: dict
        dictionary of additional arguments for the conversion function
    
    Returns
    ==================
    parameters: dict
        input parameters, but with delta_f added
    '''
        
    n = int(parameters['n'])
    total_mass = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
    parameters['delta_f'] = (203025.4467280836/total_mass)*(parameters['xi_low']**(1-1/n)*parameters['xi_high']**(1/n)-parameters['xi_low'])
    return parameters



class PriorDict(dict):
    def __init__(self, dictionary=None, filename=None, conversion_function=None, conversion_arguments=None):
        """A dictionary of priors

        Parameters
        ==========
        dictionary: Union[dict, str, None]
            If given, a dictionary to generate the prior set.
        filename: Union[str, None]
            If given, a file containing the prior to generate the prior set.
        conversion_function: func
            Function to convert between sampled parameters and constraints.
            Default is no conversion.
        """
        super(PriorDict, self).__init__()
        if isinstance(dictionary, dict):
            self.from_dictionary(dictionary)
        elif type(dictionary) is str:
            logger.debug(
                'Argument "dictionary" is a string.'
                + " Assuming it is intended as a file name."
            )
            self.from_file(dictionary)
        elif type(filename) is str:
            self.from_file(filename)
        elif dictionary is not None:
            raise ValueError("PriorDict input dictionary not understood")
        self._cached_normalizations = {}

        self.convert_floats_to_delta_functions()

        if conversion_function is not None:
            self.conversion_function = conversion_function
        else:
            self.conversion_function = self.default_conversion_function
            
        if conversion_arguments is not None:
            self.conversion_arguments = conversion_arguments
        else:
            self.conversion_arguments = None

    def evaluate_constraints(self, sample, conversion_arguments):
        if conversion_arguments is not None:
            out_sample = self.conversion_function(sample, conversion_arguments)
        else:
            out_sample = self.conversion_function(sample)
        prob = 1
        for key in self:
            if isinstance(self[key], bilby.core.prior.base.Constraint) and key in out_sample:
                prob *= self[key].prob(out_sample[key])
        return prob

    def default_conversion_function(self, sample):
        """
        Placeholder parameter conversion function.

        Parameters
        ==========
        sample: dict
            Dictionary to convert

        Returns
        =======
        sample: dict
            Same as input
        """
        return sample

    def to_file(self, outdir, label):
        """Write the prior distribution to file.

        Parameters
        ==========
        outdir: str
            output directory name
        label: str
            Output file naming scheme
        """

        check_directory_exists_and_if_not_mkdir(outdir)
        prior_file = os.path.join(outdir, "{}.prior".format(label))
        logger.debug("Writing priors to {}".format(prior_file))
        joint_dists = []
        with open(prior_file, "w") as outfile:
            for key in self.keys():
                if JointPrior in self[key].__class__.__mro__:
                    distname = "_".join(self[key].dist.names) + "_{}".format(
                        self[key].dist.distname
                    )
                    if distname not in joint_dists:
                        joint_dists.append(distname)
                        outfile.write("{} = {}\n".format(distname, self[key].dist))
                    diststr = repr(self[key].dist)
                    priorstr = repr(self[key])
                    outfile.write(
                        "{} = {}\n".format(key, priorstr.replace(diststr, distname))
                    )
                else:
                    outfile.write("{} = {}\n".format(key, self[key]))

    def _get_json_dict(self):
        self.convert_floats_to_delta_functions()
        total_dict = {key: json.loads(self[key].to_json()) for key in self}
        total_dict["__prior_dict__"] = True
        total_dict["__module__"] = self.__module__
        total_dict["__name__"] = self.__class__.__name__
        return total_dict

    def to_json(self, outdir, label):
        check_directory_exists_and_if_not_mkdir(outdir)
        prior_file = os.path.join(outdir, "{}_prior.json".format(label))
        logger.debug("Writing priors to {}".format(prior_file))
        with open(prior_file, "w") as outfile:
            json.dump(self._get_json_dict(), outfile, cls=BilbyJsonEncoder, indent=2)

    def from_file(self, filename):
        """Reads in a prior from a file specification

        Parameters
        ==========
        filename: str
            Name of the file to be read in

        Notes
        =====
        Lines beginning with '#' or empty lines will be ignored.
        Priors can be loaded from:

        - bilby.core.prior as, e.g.,    :code:`foo = Uniform(minimum=0, maximum=1)`
        - floats, e.g.,                 :code:`foo = 1`
        - bilby.gw.prior as, e.g.,      :code:`foo = bilby.gw.prior.AlignedSpin()`
        - other external modules, e.g., :code:`foo = my.module.CustomPrior(...)`

        """

        comments = ["#", "\n"]
        prior = dict()
        with ioopen(filename, "r", encoding="unicode_escape") as f:
            for line in f:
                if line[0] in comments:
                    continue
                line.replace(" ", "")
                elements = line.split("=")
                key = elements[0].replace(" ", "")
                val = "=".join(elements[1:]).strip()
                prior[key] = val
        self.from_dictionary(prior)

    @classmethod
    def _get_from_json_dict(cls, prior_dict):
        try:
            class_ = getattr(
                import_module(prior_dict["__module__"]), prior_dict["__name__"]
            )
        except ImportError:
            logger.debug(
                "Cannot import prior module {}.{}".format(
                    prior_dict["__module__"], prior_dict["__name__"]
                )
            )
            class_ = cls
        except KeyError:
            logger.debug("Cannot find module name to load")
            class_ = cls
        for key in ["__module__", "__name__", "__prior_dict__"]:
            if key in prior_dict:
                del prior_dict[key]
        obj = class_(prior_dict)
        return obj

    @classmethod
    def from_json(cls, filename):
        """Reads in a prior from a json file

        Parameters
        ==========
        filename: str
            Name of the file to be read in
        """
        with open(filename, "r") as ff:
            obj = json.load(ff, object_hook=decode_bilby_json)

        # make sure priors containing JointDists are properly handled and point
        # to the same object when required
        jointdists = {}
        for key in obj:
            if isinstance(obj[key], JointPrior):
                for name in obj[key].dist.names:
                    jointdists[name] = obj[key].dist
        # set dist for joint values so that they point to the same object
        for key in obj:
            if isinstance(obj[key], JointPrior):
                obj[key].dist = jointdists[key]

        return obj

    def from_dictionary(self, dictionary):
        mvgkwargs = {}
        for key in list(dictionary.keys()):
            val = dictionary[key]
            if isinstance(val, bilby.core.prior.base.Prior):
                continue
            elif isinstance(val, (int, float)):
                dictionary[key] = DeltaFunction(peak=val)
            elif isinstance(val, str):
                cls = val.split("(")[0]
                args = "(".join(val.split("(")[1:])[:-1]
                try:
                    dictionary[key] = DeltaFunction(peak=float(cls))
                    logger.debug("{} converted to DeltaFunction prior".format(key))
                    continue
                except ValueError:
                    pass
                if "." in cls:
                    module = ".".join(cls.split(".")[:-1])
                    cls = cls.split(".")[-1]
                else:
                    module = __name__.replace(
                        "." + os.path.basename(__file__).replace(".py", ""), ""
                    )
                try:
                    cls = getattr(import_module(module), cls, cls)
                except ModuleNotFoundError:
                    logger.error(
                        "Cannot import prior class {} for entry: {}={}".format(
                            cls, key, val
                        )
                    )
                    raise
                if key.lower() in ["conversion_function", "condition_func"]:
                    setattr(self, key, cls)
                elif isinstance(cls, str):
                    if "(" in val:
                        raise TypeError("Unable to parse prior class {}".format(cls))
                    else:
                        continue
                elif cls.__name__ in [
                    "MultivariateGaussianDist",
                    "MultivariateNormalDist",
                ]:
                    dictionary.pop(key)
                    if key not in mvgkwargs:
                        mvgkwargs[key] = cls.from_repr(args)
                elif cls.__name__ in ["MultivariateGaussian", "MultivariateNormal"]:
                    mgkwargs = {
                        item[0].strip(): cls._parse_argument_string(item[1])
                        for item in cls._split_repr(
                            ", ".join(
                                [arg for arg in args.split(",") if "dist=" not in arg]
                            )
                        ).items()
                    }
                    keymatch = re.match(r"dist=(?P<distkey>\S+),", args)
                    if keymatch is None:
                        raise ValueError(
                            "'dist' argument for MultivariateGaussian is not specified"
                        )

                    if keymatch["distkey"] not in mvgkwargs:
                        raise ValueError(
                            f"MultivariateGaussianDist {keymatch['distkey']} must be defined before {cls.__name__}"
                        )

                    mgkwargs["dist"] = mvgkwargs[keymatch["distkey"]]
                    dictionary[key] = cls(**mgkwargs)
                else:
                    try:
                        dictionary[key] = cls.from_repr(args)
                    except TypeError as e:
                        raise TypeError(
                            "Unable to parse prior, bad entry: {} "
                            "= {}. Error message {}".format(key, val, e)
                        )
            elif isinstance(val, dict):
                try:
                    _class = getattr(
                        import_module(val.get("__module__", "none")),
                        val.get("__name__", "none"),
                    )
                    dictionary[key] = _class(**val.get("kwargs", dict()))
                except ImportError:
                    logger.debug(
                        "Cannot import prior module {}.{}".format(
                            val.get("__module__", "none"), val.get("__name__", "none")
                        )
                    )
                    logger.warning(
                        "Cannot convert {} into a prior object. "
                        "Leaving as dictionary.".format(key)
                    )
                    continue
            else:
                raise TypeError(
                    "Unable to parse prior, bad entry: {} "
                    "= {} of type {}".format(key, val, type(val))
                )
        self.update(dictionary)

    def convert_floats_to_delta_functions(self):
        """Convert all float parameters to delta functions"""
        for key in self:
            if isinstance(self[key], bilby.core.prior.base.Prior):
                continue
            elif isinstance(self[key], float) or isinstance(self[key], int):
                self[key] = DeltaFunction(self[key])
                logger.debug("{} converted to delta function prior.".format(key))
            else:
                logger.debug(
                    "{} cannot be converted to delta function prior.".format(key)
                )

    def fill_priors(self, likelihood, default_priors_file=None):
        """
        Fill dictionary of priors based on required parameters of likelihood

        Any floats in prior will be converted to delta function prior. Any
        required, non-specified parameters will use the default.

        Note: if `likelihood` has `non_standard_sampling_parameter_keys`, then
        this will set-up default priors for those as well.

        Parameters
        ==========
        likelihood: bilby.likelihood.GravitationalWaveTransient instance
            Used to infer the set of parameters to fill the prior with
        default_priors_file: str, optional
            If given, a file containing the default priors.


        Returns
        =======
        prior: dict
            The filled prior dictionary

        """

        self.convert_floats_to_delta_functions()

        missing_keys = set(likelihood.parameters) - set(self.keys())

        for missing_key in missing_keys:
            if not self.test_redundancy(missing_key):
                default_prior = create_default_prior(missing_key, default_priors_file)
                if default_prior is None:
                    set_val = likelihood.parameters[missing_key]
                    logger.warning(
                        "Parameter {} has no default prior and is set to {}, this"
                        " will not be sampled and may cause an error.".format(
                            missing_key, set_val
                        )
                    )
                else:
                    self[missing_key] = default_prior

        for key in self:
            self.test_redundancy(key)

    def sample(self, size=None):
        """Draw samples from the prior set

        Parameters
        ==========
        size: int or tuple of ints, optional
            See numpy.random.uniform docs

        Returns
        =======
        dict: Dictionary of the samples
        """
        return self.sample_subset_constrained(keys=list(self.keys()), size=size)

    def sample_subset_constrained_as_array(self, keys=iter([]), size=None):
        """Return an array of samples

        Parameters
        ==========
        keys: list
            A list of keys to sample in
        size: int
            The number of samples to draw

        Returns
        =======
        array: array_like
            An array of shape (len(key), size) of the samples (ordered by keys)
        """
        samples_dict = self.sample_subset_constrained(keys=keys, size=size)
        samples_dict = {key: np.atleast_1d(val) for key, val in samples_dict.items()}
        samples_list = [samples_dict[key] for key in keys]
        return np.array(samples_list)

    def sample_subset(self, keys=iter([]), size=None):
        """Draw samples from the prior set for parameters which are not a DeltaFunction

        Parameters
        ==========
        keys: list
            List of prior keys to draw samples from
        size: int or tuple of ints, optional
            See numpy.random.uniform docs

        Returns
        =======
        dict: Dictionary of the drawn samples
        """
        self.convert_floats_to_delta_functions()
        samples = dict()
        for key in keys:
            if isinstance(self[key], bilby.core.prior.base.Constraint):
                continue
            elif isinstance(self[key], bilby.core.prior.base.Prior):
                samples[key] = self[key].sample(size=size)
            else:
                logger.debug("{} not a known prior.".format(key))
        return samples

    @property
    def non_fixed_keys(self):
        keys = self.keys()
        keys = [k for k in keys if isinstance(self[k], bilby.core.prior.base.Prior)]
        keys = [k for k in keys if self[k].is_fixed is False]
        keys = [k for k in keys if k not in self.constraint_keys]
        return keys

    @property
    def fixed_keys(self):
        return [
            k for k, p in self.items() if (p.is_fixed and k not in self.constraint_keys)
        ]

    @property
    def constraint_keys(self):
        return [k for k, p in self.items() if isinstance(p, bilby.core.prior.base.Constraint)]

    def sample_subset_constrained(self, keys=iter([]), size=None, conversion_arguments=None):
        if size is None or size == 1:
            while True:
                sample = self.sample_subset(keys=keys, size=size)
                if self.evaluate_constraints(sample,self.conversion_arguments):
                    return sample
        else:
            needed = np.prod(size)
            for key in keys.copy():
                if isinstance(self[key], bilby.core.prior.base.Constraint):
                    del keys[keys.index(key)]
            all_samples = {key: np.array([]) for key in keys}
            _first_key = list(all_samples.keys())[0]
            while len(all_samples[_first_key]) < needed:
                samples = self.sample_subset(keys=keys, size=needed)
                keep = np.array(self.evaluate_constraints(samples), dtype=bool)
                for key in keys:
                    all_samples[key] = np.hstack(
                        [all_samples[key], samples[key][keep].flatten()]
                    )
            all_samples = {
                key: np.reshape(all_samples[key][:needed], size) for key in keys
            }
            return all_samples

    def normalize_constraint_factor(
        self, keys, min_accept=10000, sampling_chunk=50000, nrepeats=10
    ):
        if keys in self._cached_normalizations.keys():
            return self._cached_normalizations[keys]
        else:
            factor_estimates = [
                self._estimate_normalization(keys, min_accept, sampling_chunk)
                for _ in range(nrepeats)
            ]
            factor = np.mean(factor_estimates)
            if np.std(factor_estimates) > 0:
                decimals = int(-np.floor(np.log10(3 * np.std(factor_estimates))))
                factor_rounded = np.round(factor, decimals)
            else:
                factor_rounded = factor
            self._cached_normalizations[keys] = factor_rounded
            return factor_rounded

    def _estimate_normalization(self, keys, min_accept, sampling_chunk):
        samples = self.sample_subset(keys=keys, size=sampling_chunk)
        keep = np.atleast_1d(self.evaluate_constraints(samples))
        if len(keep) == 1:
            self._cached_normalizations[keys] = 1
            return 1
        all_samples = {key: np.array([]) for key in keys}
        while np.count_nonzero(keep) < min_accept:
            samples = self.sample_subset(keys=keys, size=sampling_chunk)
            for key in samples:
                all_samples[key] = np.hstack([all_samples[key], samples[key].flatten()])
            keep = np.array(self.evaluate_constraints(all_samples), dtype=bool)
        factor = len(keep) / np.count_nonzero(keep)
        return factor

    def prob(self, sample, **kwargs):
        """

        Parameters
        ==========
        sample: dict
            Dictionary of the samples of which we want to have the probability of
        kwargs:
            The keyword arguments are passed directly to `np.prod`

        Returns
        =======
        float: Joint probability of all individual sample probabilities

        """
        prob = np.prod([self[key].prob(sample[key]) for key in sample], **kwargs)

        return self.check_prob(sample, prob)

    def check_prob(self, sample, prob):
        ratio = self.normalize_constraint_factor(tuple(sample.keys()))
        if np.all(prob == 0.0):
            return prob * ratio
        else:
            if isinstance(prob, float):
                if self.evaluate_constraints(sample):
                    return prob * ratio
                else:
                    return 0.0
            else:
                constrained_prob = np.zeros_like(prob)
                keep = np.array(self.evaluate_constraints(sample), dtype=bool)
                constrained_prob[keep] = prob[keep] * ratio
                return constrained_prob

    def ln_prob(self, sample, axis=None, normalized=True):
        """

        Parameters
        ==========
        sample: dict
            Dictionary of the samples of which to calculate the log probability
        axis: None or int
            Axis along which the summation is performed
        normalized: bool
            When False, disables calculation of constraint normalization factor
            during prior probability computation. Default value is True.

        Returns
        =======
        float or ndarray:
            Joint log probability of all the individual sample probabilities

        """
        ln_prob = np.sum([self[key].ln_prob(sample[key]) for key in sample], axis=axis)
        return self.check_ln_prob(sample, ln_prob,
                                  normalized=normalized)

    def check_ln_prob(self, sample, ln_prob, normalized=True):
        if normalized:
            ratio = self.normalize_constraint_factor(tuple(sample.keys()))
        else:
            ratio = 1
        if np.all(np.isinf(ln_prob)):
            return ln_prob
        else:
            if isinstance(ln_prob, float):
                if self.evaluate_constraints(sample):
                    return ln_prob + np.log(ratio)
                else:
                    return -np.inf
            else:
                constrained_ln_prob = -np.inf * np.ones_like(ln_prob)
                keep = np.array(self.evaluate_constraints(sample), dtype=bool)
                constrained_ln_prob[keep] = ln_prob[keep] + np.log(ratio)
                return constrained_ln_prob

    def cdf(self, sample):
        """Evaluate the cumulative distribution function at the provided points

        Parameters
        ----------
        sample: dict, pandas.DataFrame
            Dictionary of the samples of which to calculate the CDF

        Returns
        -------
        dict, pandas.DataFrame: Dictionary containing the CDF values

        """
        return sample.__class__(
            {key: self[key].cdf(sample) for key, sample in sample.items()}
        )

    def rescale(self, keys, theta):
        """Rescale samples from unit cube to prior

        Parameters
        ==========
        keys: list
            List of prior keys to be rescaled
        theta: list
            List of randomly drawn values on a unit cube associated with the prior keys

        Returns
        =======
        list: List of floats containing the rescaled sample
        """
        from matplotlib.cbook import flatten

        return list(
            flatten([self[key].rescale(sample) for key, sample in zip(keys, theta)])
        )

    def test_redundancy(self, key, disable_logging=False):
        """Empty redundancy test, should be overwritten in subclasses"""
        return False

    def test_has_redundant_keys(self):
        """
        Test whether there are redundant keys in self.

        Returns
        =======
        bool: Whether there are redundancies or not
        """
        redundant = False
        for key in self:
            if isinstance(self[key], bilby.core.prior.base.Constraint):
                continue
            temp = self.copy()
            del temp[key]
            if temp.test_redundancy(key, disable_logging=True):
                logger.warning(
                    "{} is a redundant key in this {}.".format(
                        key, self.__class__.__name__
                    )
                )
                redundant = True
        return redundant

    def copy(self):
        """
        We have to overwrite the copy method as it fails due to the presence of
        defaults.
        """
        return self.__class__(dictionary=dict(self))
