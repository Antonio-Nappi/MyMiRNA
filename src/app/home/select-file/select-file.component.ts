import { Component, OnInit, Input } from '@angular/core';
import { ModalController } from '@ionic/angular';
import { NgForm } from '@angular/forms';

@Component({
  selector: 'app-select-file',
  templateUrl: './select-file.component.html',
  styleUrls: ['./select-file.component.scss'],
})

export class SelectFileComponent implements OnInit {

  @Input() files = new Array();

  public times = new Array();

  constructor(private modalCtrl: ModalController) { }

  ngOnInit() {
  }
/*
  onSubmit(form: NgForm, num: number) {
    if (form.valid) {
      for (let i = 0; i < num; i++) {
        this.times.push('file' + i);
      }
      document.getElementById('input').setAttribute('disabled', 'true');
      document.getElementById('button').setAttribute('disabled', 'true');
      document.getElementById('confirm').setAttribute('disabled', 'false');
    }
  }

  onSubmitModal(form: NgForm) {
    for (let i=0; i < this.times.length; i++) {
      this.files.push(form.value[this.times.pop()]);
    }
    this.modalCtrl.dismiss(this.files, 'confirm');
  }

  onCloseModal() {
    this.modalCtrl.dismiss(null, 'close');
  }

*/

}
