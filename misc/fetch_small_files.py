import sys
import requests

def fetch_file(file_url):
    res = requests.get(file_url)
    if res.status_code != 200:
        sys.stderr.write("Failed to fetch {}".format(file_url))
    file_name = res.headers.get("content-disposition")
    file_content = res.content
    return file_name, file_content

def write_file(file_name, file_content):
    with open(file_name, "wb") as file_fd:
        fd.write(file_content)

def parse_urls(urls_file_name):
    with open(urls_file_name, "r") as urls_fd:
        urls = urls_fd.readlines()
    return urls

def main(args):
    urls = parse_urls(args[0])
    for url in urls:
        file_name, file_content = fetch_file(url)
        write_file(file_name, file_content)

if __name__ == "__main__":
    main(sys.argv[1:])