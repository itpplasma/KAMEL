"""Public exception hierarchy for the KIM Python interface."""


class KimError(Exception):
    """Base class for errors raised by the KIM Python interface."""


class ConfigurationError(KimError):
    """Raised when a KIM configuration cannot be constructed or loaded."""


class ProfileError(KimError):
    """Raised when radial profile input is missing, malformed, or inconsistent."""
