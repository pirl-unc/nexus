isON_pipeline.sh is a modified version of the [isON_pipeline.sh](https://github.com/aljpetri/isONform/blob/master/isON_pipeline.sh). It modified the script in the following ways:

### 1. Added `die` function
```sh
function die {
    echo "$@" >&2
    exit 1
}
```
Every "Missing parameter" branch called die, but it was never defined anywhere. 
Under set -e those branches aborted with die: command not found and exit 127, 
swallowing the message they were meant to print. 
Now all four report the actual missing parameter and exit 1.


### 2. Added `--iso_abundance` validation
```sh
elif [[ -z $iso_abundance ]]; then
    usage
    die "Missing parameter --iso_abundance"
```


### 3. `python isONform_parallel --h` is now `isONform_parallel --h`
isONform_parallel is an installed console script on PATH, so python isONform_parallel 
makes Python look for a file by that name in the working directory and fail. 