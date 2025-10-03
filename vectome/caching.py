from platformdirs import user_cache_dir

from . import app_name, __author__, __version__

CACHE_DIR = user_cache_dir(
    app_name, 
    __author__,
)
