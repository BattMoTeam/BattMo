(with-current-buffer (find-file-noselect "temp.sh")
  (erase-buffer)
  (cl-loop for filename in (directory-files "." nil "pdf$") do
           (let ((pngfilename (string-replace "pdf" "png" filename)))
             (insert "convert -density 1000 " filename " " pngfilename "\n")
             (insert "mv -f " pngfilename " ../\n")))
  (save-buffer)
  (shell-command "sh temp.sh")
  )
         

