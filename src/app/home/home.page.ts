import { Component, OnInit } from '@angular/core';
import { NgForm } from '@angular/forms';
import { RequestService } from './requests.service';
import { ModalController } from '@ionic/angular';
import { SelectFileComponent } from './select-file/select-file.component';

@Component({
  selector: 'app-home',
  templateUrl: 'home.page.html',
  styleUrls: ['home.page.scss'],
})
export class HomePage implements OnInit {

  private paths = new Array();
  private multiqc: any;

  constructor(private request: RequestService, private modalCtrl: ModalController) { }

  ngOnInit() {
    document.getElementById('card2').setAttribute('disabled', 'true');
    document.getElementById('card3').setAttribute('disabled', 'true');
    document.getElementById('card4').setAttribute('disabled', 'true');
    document.getElementById('card5').setAttribute('disabled', 'true');
    document.getElementById('card6').setAttribute('disabled', 'true');
    document.getElementById('card7').setAttribute('disabled', 'true');
  }

  onSubmit(form: NgForm, name: string) {
    if (form.valid) {
      if (name === 'f1') {
        const body = JSON.stringify({
          qual: form.value['qual'],
          adapter: form.value['adapter']
        });
        console.log(body);
        this.multiqc = this.request.trimmingStep(body);
        document.getElementById('card2').setAttribute('disabled', 'false');
      } else if (name === 'f2') {
        const body = JSON.stringify({
          multimap: form.value['tresh'],
          cores: form.value['core'],
          p_value: form.value['p_value'],
          p_value_adjusted: form.value['p_value_adjusted'],
          log_2_fold: form.value['log']
        });
        console.log(body);
        document.getElementById('card3').setAttribute('disabled', 'false');
      } else if (form.name === 'f3') {
        const body = JSON.stringify({
          multimap: form.value['mismatches'],
          cores: form.value['core'],
          p_value: form.value['p_value'],
          p_value_adjusted: form.value['p_value_adjusted'],
          log_2_fold: form.value['log']
        });
        console.log(body);
        document.getElementById('card4').setAttribute('disabled', 'false');
      } else if (name === 'f4') {
        const body = JSON.stringify({
          multimap: form.value['mismatches'],
          cores: form.value['core'],
          p_value: form.value['p_value'],
          p_value_adjusted: form.value['p_value_adjusted'],
          log_2_fold: form.value['log']
        });
        console.log(body);
        document.getElementById('card5').setAttribute('disabled', 'false');
      } else if (name === 'f5') {
        const body = JSON.stringify({
          multimap: form.value['mismatches'],
          cores: form.value['core'],
          p_value: form.value['p_value'],
          p_value_adjusted: form.value['p_value_adjusted'],
          log_2_fold: form.value['log']
        });
        console.log(body);
        document.getElementById('card6').setAttribute('disabled', 'false');
      } else if (name === 'f6') {
        const body = JSON.stringify({
          multimap: form.value['mismatches'],
          cores: form.value['core'],
          p_value: form.value['p_value'],
          p_value_adjusted: form.value['p_value_adjusted'],
          log_2_fold: form.value['log']
        });
        console.log(body);
        document.getElementById('card7').setAttribute('disabled', 'false');
      } else {
        const body = JSON.stringify({
          multimap: form.value['mismatches'],
          cores: form.value['core'],
          p_value: form.value['p_value'],
          p_value_adjusted: form.value['p_value_adjusted'],
          log_2_fold: form.value['log']
        });
        console.log(body);
      }
    }
  }
  /*
    openModal() {
      this.modalCtrl.create({
        component: SelectFileComponent,
        componentProps: { files: this.paths }
      }).then(modalEl => {
        modalEl.present();
        return modalEl.onDidDismiss();
      }).then(resData => {
        this.paths = resData.data;
      });
    }
    */

}
