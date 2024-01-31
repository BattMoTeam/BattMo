;;; Utils for building BattMo doc

;; Use user-login-name to get user name
;; (cond '(compare-strings user-login-name "xavier")

(pcase (user-login-name)
  ("xavier" (progn
              (pyvenv-activate "~/Python/battmodoc-env/")
              (setq docdir "/home/xavier/Matlab/Projects/battmo/Documentation/")
              (setq rootdir "/home/xavier/Matlab/Projects/battmo/")
              (setq utilsdir "/home/xavier/Matlab/Projects/battmo/Documentation/utils/")
              (setq testdir "/home/xavier/Matlab/Projects/battmo-doc-test/")))
  ("august" (progn
              (setq docdir "/home/august/Projects/Battery/BattMo/Documentation/")
              (setq utilsdir "/home/august/Projects/Battery/BattMo/Documentation/utils/")
              (setq testdir "/home/august/Projects/Battery/battmo-doc-test/")

              )
   )
  )

(defun battmodoc-local-open ()
  "Open locally built documentation in browser"
  (interactive)
  (browse-url (concat docdir "_build/html/index.html"))
  )

(defun battmodoc-build ()
  "Build BattMo documentation"
  (interactive)
  (let* ((default-directory docdir)
         (outputbuffer (get-buffer-create "*buildoutput*"))
         )
    (pop-to-buffer outputbuffer)
    (erase-buffer)
    (start-process "battmo-build" outputbuffer "make" "html")
    )
  )

(defun battmodoc-build-examples ()
  "Run python publish script"
  (interactive)
  (let ((outputbuffer (get-buffer-create "*publishoutput*"))
        (directory utilsdir))
    (pop-to-buffer outputbuffer)
    (erase-buffer)
    (start-process "battmo-publish" outputbuffer "python" (concat directory "buildPublishedExamples.py"))
    )
  )

(defun battmodoc-publish ()
  "Copy and paste build directory in doc repo and push the result. Then open in a browser the action page in the repo and the published page"
  (interactive)
  (let ((default-directory docdir))
    (shell-command (concat "cp -rf _build/html/* " testdir))
    )
  (let ((default-directory testdir))
    (shell-command "git add *" )
    (shell-command "git commit -m \"update in doc\"")
    (shell-command "git push -f" )
    )
  (browse-url "https://github.com/BattMoTeam/battmo-doc-test/actions")
  (browse-url "https://battmoteam.github.io/battmo-doc-test/")
  )

(defun get-broken-links (url)
  
  (let* ((base-url "https://battmoteam.github.io/BattMo/")
         (url (concat base-url url))
         (wget-buffer (pop-to-buffer "*battmo check link*"))
         (broken-links nil))

    ;;  The following are the basic flags you'll need for the wget commandd:
    ;;
    ;;  --spider stops wget from downloading the page.
    ;;  -r makes wget recursively follow each link on the page.
    ;;  -nd, short for --no-directories, prevents wget from creating a hierarchy of directories on your server (even when it is configured to spider only).
    ;;  -nv, short for --no-verbose, stops wget from outputting extra information that is unnecessary for identifying broken links.
    ;;
    ;;  The following are optional parameters which you can use to customize your search:
    ;;  
    ;;  -H, short for --span-hosts,makes wget crawl to subdomains and domains other than the primary one (i.e. external sites).
    ;;  -l 1 is short for --level. By default, wget crawls up to five levels deep from the initial URL, but here we set it to one. You may need to play with this parameter depending on the organization of your website.
    ;;  -w 2, short for --wait, instructs wget to wait 2 seconds between requests to avoid bombarding the server, minimizing any performance impact.
    ;;  -o run1.log saves wget’s output to a file called run1.log instead of displaying it in your terminal.
    (with-current-buffer wget-buffer
      (erase-buffer)
      (call-process  "wget" nil wget-buffer t "--spider" "-r" "-nd" "-nv" "-H" "-l 1" "-w 0.01" url)
      (beginning-of-buffer)
      (while (re-search-forward (rx "Remote file does not exist -- broken link!!!") nil t)
        (save-excursion
          (previous-line)
          (push (buffer-substring-no-properties (line-beginning-position) (line-end-position)) broken-links)
          )
        )
      )
    broken-links
    )
  )

(defun collect-all-broken-links ()
  (let* ((build-dir (concat docdir "_build/html/"))
         (html-files (cl-loop for file in (directory-files-recursively build-dir (rx ".html" eol)) collect (string-replace build-dir "" file)))
         )
    (cl-loop for file in html-files collect (list file (get-broken-links file)))
    ))

(setq results (collect-all-broken-links))

(with-current-buffer (find-file-noselect "temp.org")
  (erase-buffer)
  (cl-loop for res in results do
           (pcase res
             (`(,file ,links) (if links
                                  (progn
                                    (insert "* " file "\n")
                                    (cl-loop for link in links do (insert " - " link "\n"))
                                    )
                                )
              )
             )
           )
  (save-buffer)
  )

;; (setq res (check-broken-links "basicusage.html"))

(defun convert-to-attribute ()
  (interactive)
  (re-search-forward (rx (group (1+ (any word))) (0+ space) "%" (0+ space) (group (0+ nonl))))
  (let ((att (match-string 1))
        (com (match-string 2))
        )
    (beginning-of-line)
    (kill-line)
    (insert "   .. attribute:: " att "\n\n")
    (insert "      " com "\n")
    )
  )

(defun battmodoc-matlab-include-link ()
  (interactive)
  (let* ((url "https://github.com/BattMoTeam/BattMo")
        (branch-name "dev")
        (fullfilename (battmo-counsel-fzf "" rootdir))
        (filename (file-name-nondirectory fullfilename))
        (p (make-marker)))
    (insert (concat "<" url "/blob/" branch-name "/" fullfilename "  "))
    (set-marker p (point))
    (insert (concat filename ">"))
    (goto-char p)
    )
  )

(defun battmo-counsel-fzf (initial-input initial-directory)
  (counsel-require-program counsel-fzf-cmd)
  (setq counsel--fzf-dir  initial-directory)
  (ivy-read "file to link:"
            #'counsel-fzf-function
            :initial-input initial-input
            :re-builder #'ivy--regex-fuzzy
            :dynamic-collection t
            :action #'identity
            :caller 'counsel-fzf))


